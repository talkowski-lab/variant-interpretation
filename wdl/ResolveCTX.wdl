version 1.0

# IMPORT
import "Structs.wdl"
import "TinyResolve.wdl"
import "PEevidence.wdl" as PEevidence

# WORKFLOW DEFINITION
workflow ResolveCTX{
    input{
        Array[String] batches

        File vcf
        String prefix
        String docker_path
        Array[String] samples
        File manta_vcf_tar
        File cytoband
        Array[File] discfile
        File mei_bed
        Int samples_per_shard = 25
        String sv_pipeline_docker
        String linux_docker

        File batches_pe
        File sample_batch
        String docker_pe_evidence

        RuntimeAttr? runtime_attr_resolve
        RuntimeAttr? runtime_attr_untar
        RuntimeAttr? runtime_attr_reformat
        RuntimeAttr? runtime_attr_override_merge
        RuntimeAttr? runtime_attr_subset_ctx_vcf
        RuntimeAttr? runtime_attr_override_get_raw_only
        RuntimeAttr? runtime_attr_override_ref_raw_pe
        RuntimeAttr? runtime_attr_subset_sample_roi
        RuntimeAttr? runtime_attr_subset_pe_evidence
    }

    ##Get CTXs from main GATK-SV output
    call CtxVcf2Bed{
        input:
            vcf = vcf,
            prefix = prefix,
            docker_path = docker_path,
            runtime_attr_override = runtime_attr_subset_ctx_vcf
    }

    ##Get CTXx from Manta raw evidence after running TinyResolve (per batch)
    scatter (batch in batches){

        call TinyResolve.TinyResolve as TinyResolve{
            input:
                samples = samples,
                manta_vcf_tar = manta_vcf_tar,
                cytoband = cytoband,
                discfile = discfile,
                mei_bed = mei_bed,
                samples_per_shard = samples_per_shard,
                sv_pipeline_docker = sv_pipeline_docker,
                linux_docker = linux_docker,
                runtime_attr_resolve = runtime_attr_resolve,
                runtime_attr_untar = runtime_attr_untar
        }

        scatter (file in TinyResolve.tloc_manta_vcf){
            String file_idx = basename(file, ".tbi")

            call reformatTinyResolve{
                input:
                    input_vcf=file,
    #                input_vcf_idx=file_idx,
                    docker = docker_path,
                    runtime_attr_override = runtime_attr_reformat
            }
        }

        call mergeTinyResolve{
            input:
                input_beds = reformatTinyResolve.ctx_raw_for_vcf_ovl,
                docker = docker_path,
                prefix = prefix,
                runtime_attr_override = runtime_attr_override_merge
        }

        call getRawOnlyCTX{
            input:
                ctx_vcf_input = CtxVcf2Bed.ctx_vcf_for_raw_ovl,
                ctx_raw_input = mergeTinyResolve.merged_ctx_raw_for_vcf_ovl,
                docker = docker_path,
                prefix = prefix,
                runtime_attr_override = runtime_attr_override_get_raw_only
        }

        call reformatRawOnlyForPE{
            input:
                ctx_raw_input = getRawOnlyCTX.raw_only,
                docker = docker_path,
                prefix = prefix,
                runtime_attr_override = runtime_attr_override_ref_raw_pe
        }
    }

    call mergeVcfRawForPE{
        input:
            input_raw = reformatRawOnlyForPE.raw_only_for_pe,
            input_vcf = CtxVcf2Bed.ctx_vcf_for_pe,
            docker = docker_path,
            prefix = prefix,
            runtime_attr_override = runtime_attr_override_merge
    }

    ###TO DO AFTER MERGING ACROSS BATCHES???????
    File samples_for_pe = mergeVcfRawForPE.samples_with_ctx_for_pe
    Array[String] samples_pe = transpose(read_tsv(samples_for_pe))[0]

    scatter (sample in samples_pe){

        call PEevidence.subset_sample_roi as subset_sample_roi{
            input:
                sample = sample,
                regions = mergeVcfRawForPE.merged_vcf_raw_for_pe,
                docker_path = docker_pe_evidence,
                runtime_attr_override = runtime_attr_subset_sample_roi
        }

        call PEevidence.subset_pe_evidence as subset_pe_evidence{
            input:
                sample_bed = subset_sample_roi.sample_roi,
                sample_batch = sample_batch,
                batches_pe = batches_pe,
                sample = sample,
                docker_path = docker_pe_evidence,
                runtime_attr_override = runtime_attr_subset_pe_evidence
        }
    }

    call mergeSamplePE{
        input:
            input_beds = subset_pe_evidence.sample_pe,
            docker = docker_path,
            prefix = prefix,
            runtime_attr_override = runtime_attr_override_merge
    }

    output{
        File ctx_vcf_for_pe = CtxVcf2Bed.ctx_vcf_for_pe
        File ctx_vcf_for_raw_ovl = CtxVcf2Bed.ctx_vcf_for_raw_ovl
        File tinyresolved_merged = mergeTinyResolve.merged_ctx_raw_for_vcf_ovl
        File raw_vcf_ovl_bed = getRawOnlyCTX.raw_vcf_ovl
        File raw_only_manta_bed = getRawOnlyCTX.raw_only
        File ctx_raw_for_pe = reformatRawOnlyForPE.raw_only_for_pe
#        File ctx_vcf_raw_merged_for_pe = mergeVcfRawForPE.merged_vcf_raw_for_pe
        File merged_pe_evidence = mergeSamplePE.merged_PE_evidence
    }
}

# TASK DEFINITIONS
task CtxVcf2Bed {
  input {
    File vcf
    String prefix
    String docker_path
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output{
    File ctx_vcf_for_pe = "~{prefix}.ctx.refForPE.split.bed.gz"
    File ctx_vcf_for_raw_ovl = "~{prefix}.ctx.refForOverlap.split.bed.gz"
  }

  command <<<
    set -euo pipefail
    bcftools view ~{vcf} | grep -E "^#|SVTYPE=CTX" | bgzip -c > ~{prefix}.ctx.vcf.gz
    svtk vcf2bed -i ALL --include-filters ~{prefix}.ctx.vcf.gz - | bgzip -c > ~{prefix}.ctx.bed.gz

    zcat ~{prefix}.ctx.bed.gz | awk 'BEGIN { OFS="\t" } {print $1, $2, $8, $12, $6, $36, $4}'| tail -n+2 | bgzip -c > ~{prefix}.ctx.refForPE.bed.gz
    zcat ~{prefix}.ctx.refForPE.bed | awk 'BEGIN { OFS="\t" } {split($5,a,",");for(i in a)if(!seen[a[i]]++)print $1, $2, $3, $4, a[i], $6, $7}' | bgzip -c > ~{prefix}.ctx.refForPE.split.bed.gz
    zcat ~{prefix}.ctx.refForPE.split.bed | awk 'BEGIN { OFS="\t" } {print $1"_"$5, $2-20, $2+20, substr($0, index($0,$3))}' | bgzip -c > ~{prefix}.ctx.refForOverlap.split.bed.gz

  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_path
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task reformatTinyResolve{
    input{
        File input_vcf
#        File input_vcf_idx
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(input_vcf, ".vcf.gz")

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible_tries: 3,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File ctx_raw_for_vcf_ovl = "~{prefix}_raw.ctx.refForOvl.bed.gz"
    }
    command <<<
        set -euo pipefail

        svtk vcf2bed -i ALL --include-filters ~{input_vcf} ~{prefix}.bed
        awk 'BEGIN { OFS="\t" } /^#/; $5 == "CTX"{print $1"_"$6,$2-20,$3+20,substr($0, index($0,$4))}' ~{prefix}.bed | \
            bgzip -c > ~{prefix}_raw.ctx.refForOvl.bed.gz
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}

task mergeTinyResolve{
    input{
        Array [File] input_beds
        String docker
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible_tries: 3,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File merged_ctx_raw_for_vcf_ovl = "~{prefix}_raw_ctx.bed.gz"
    }
    command <<<
        set -euo pipefail
        zcat ~{sep=' ' input_beds} | \
            grep -v "^#chrom" | \
            sort -k 1,1 -k2,2n | \
            bgzip -c > ~{prefix}_raw_ctx.bed.gz
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}

task getRawOnlyCTX{
    input{
        File ctx_vcf_input
        File ctx_raw_input
        String docker
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible_tries: 3,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File raw_only = "~{prefix}_raw_only.bed.gz"
        File raw_vcf_ovl = "~{prefix}_raw_vcf_ovl.bed.gz"
    }
    command <<<
        set -euo pipefail
        bedtools intersect -a ~{ctx_raw_input} -b ~{ctx_vcf_input} |
            bgzip -c > ~{prefix}_raw_vcf_ovl.bed.gz
        bedtools intersect -v -a ~{ctx_raw_input} -b ~{ctx_vcf_input} |
            bgzip -c > ~{prefix}_raw_only.bed.gz
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}

task reformatRawOnlyForPE{
    input{
        File ctx_raw_input
        String docker
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible_tries: 3,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File raw_only_for_pe = "~{prefix}_raw_only.refForPE.bed.gz"
    }
    command <<<
        set -euo pipefail

        Rscript /src/variant-interpretation/scripts/reformatRawForPE.R ~{ctx_raw_input} ~{prefix}_raw_only.refForPE.bed
        bgzip ~{prefix}_raw_only.refForPE.bed
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}

task mergeVcfRawForPE{
    input{
        File input_raw
        File input_vcf
        String docker
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible_tries: 3,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File merged_vcf_raw_for_pe = "~{prefix}_vcf_raw_forPE.bed.gz"
        File samples_with_ctx_for_pe = "~{prefix}_samples_forPE.txt"
    }
    command <<<
        set -euo pipefail
        zcat ~{input_raw} ~{input_vcf}| \
            sort -k 1,1 -k2,2n | \
            bgzip -c > ~{prefix}_vcf_raw_forPE.bed.gz

        zcat ~{prefix}_vcf_raw_forPE.bed.gz | cut -f5 | sort -u > ~{prefix}_samples_forPE.txt
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}

task mergeSamplePE{
    input{
        Array[File] input_beds
        String docker
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible_tries: 3,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File merged_PE_evidence = "~{prefix}_PEevidence.bed.gz"
    }
    command <<<
        set -euo pipefail
        zcat ~{sep=' ' input_beds} |
            bgzip -c > ~{prefix}_PEevidence.bed.gz
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}