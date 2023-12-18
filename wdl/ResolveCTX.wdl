version 1.0

# IMPORT
import "Structs.wdl"
import "TinyResolve.wdl"
import "PEevidence.wdl"

# WORKFLOW DEFINITION
workflow ResolveCTX{
    input{
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
        RuntimeAttr? runtime_attr_resolve
        RuntimeAttr? runtime_attr_untar
        RuntimeAttr? runtime_attr_reformat
#        RuntimeAttr? runtime_attr_override_merge

        RuntimeAttr? runtime_attr_subset_ctx_vcf
    }

    call CtxVcf2Bed{
        input:
            vcf = vcf,
            prefix = prefix,
            docker_path = docker_path,
            runtime_attr_override = runtime_attr_subset_ctx_vcf
    }

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
                input_vcf_idx=file_idx,
                docker = docker_path,
                runtime_attr_override = runtime_attr_reformat
        }
    }

    call mergeTinyResolve{
        input:
            input_beds=reformatTinyResolve.ref_tiny_resolve,
            docker = docker_path,
            prefix = prefix,
            runtime_attr_override = runtime_attr_override_merge
    }


#    call ReformatTinyResolve
#    call MergeVCFtinyResolve
#    call PEevidence

    output{
        File ctx_ref_bed = CtxVcf2Bed.ctx_bed
        File tinyresolved_merged = mergeTinyResolve.merged_manta
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
    File ctx_bed = "~{prefix}.ctx.ref.split.bed"
  }

  command <<<
    set -euo pipefail
    bcftools view ~{vcf} | grep -E "^#|SVTYPE=CTX" | bgzip -c > ~{prefix}.ctx.vcf.gz
    svtk vcf2bed -i ALL --include-filters ~{prefix}.ctx.vcf.gz - | bgzip -c > ~{prefix}.ctx.bed.gz

    zcat ~{prefix}.ctx.bed.gz | awk '{print $1"\t"$2"\t"$8"\t"$12"\t"$6"\t"$36"\t"$4}'| tail -n+2 > ~{prefix}.ctx.ref.bed

    awk '{split($5,a,",");for(i in a)if(!seen[a[i]]++)print $1"\t"$2"\t"$3"\t"$4"\t"a[i]"\t"$6"\t"$7}' ~{prefix}.ctx.ref.bed > ~{prefix}.ctx.ref.split.bed
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
        File input_vcf_idx
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
        File ref_tiny_resolve = "~{prefix}_tloc.bed.gz"
    }
    command <<<
        set -euo pipefail

        svtk vcf2bed -i ALL --include-filters ~{input_vcf} ~{prefix}.bed
        awk '$5 == "CTX"{print $0}' ~{prefix}_tloc.bed | bgzip -c > ~{prefix}_tloc.bed
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
        File merged_manta = "tinyresolve_merged_tlocs.bed.gz"
    }
    command <<<
        set -euo pipefail
        zcat ~{sep=' ' input_beds} | bgzip -c > tinyresolve_merged_tlocs.bed.gz
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