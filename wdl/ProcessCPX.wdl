version 1.0

# IMPORT
## modified from ResolveCTX.wdl

import "Structs.wdl"
import "TinyResolveCPX.wdl"
import "PEevidence.wdl" as PEevidence

# WORKFLOW DEFINITION
workflow ResolveCTX{
    input{
        String prefix # batchid
        String docker_path
        Array[String] samples
        File manta_vcf_tar # batch std manta tarball
        File cytoband
        File cytoband_idx
        Array[File] discfile # per sample pesr_disc
        File mei_bed
        Int samples_per_shard = 25
        String sv_pipeline_docker
        String linux_docker

        File batches_pe # key file with batch ID + batch level pe.txt.gz + pe.txt.gz.tbi
        File sample_batch # key file with batch ID + sample ID
        String docker_pe_evidence

        RuntimeAttr? runtime_attr_resolve
        RuntimeAttr? runtime_attr_untar
        RuntimeAttr? runtime_attr_extract
    }

    call TinyResolveCPX.TinyResolveCPX as TinyResolveCPX{
        input:
            samples = samples,
            manta_vcf_tar = manta_vcf_tar,
            cytoband = cytoband,
            cytoband_idx = cytoband_idx
            discfile = discfile,
            mei_bed = mei_bed,
            samples_per_shard = samples_per_shard,
            sv_pipeline_docker = sv_pipeline_docker,
            linux_docker = linux_docker,
            runtime_attr_resolve = runtime_attr_resolve,
            runtime_attr_untar = runtime_attr_untar
    }

    scatter (file in TinyResolveCPX.cpx_manta_unresolved_vcf){
        String file_idx = basename(file, ".tbi")

        call extract_complex{
            input:
                input_vcf=file,
#                input_vcf_idx=file_idx,
                docker = docker_path,
                runtime_attr_override = runtime_attr_extract
        }
    }

    call clusterCPX{
        input:
            input_beds = extract_complex.cpx_formatted,
            docker = docker_path,
            prefix = prefix,
            runtime_attr_override = runtime_attr_override_cluster
    }

    output{
        File cluster_bed = clusterCPX.svtk_bedcluster

    }
}

# TASK DEFINITIONS

task extract_complex{
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
        File cpx_formatted = "~{input_vcf}_complex_events_formatted.bed"
    }
    command <<<
        set -euo pipefail

        # convert to bed file
        svtk vcf2bed -i ALL --include-filters ~{input_vcf} - | bgzip -c > ~{input_vcf}.bed.gz

        #extract multiple events
        zcat ~{input_vcf}.bed.gz | awk '{print $18}' | sort | uniq -c | \
        awk '{if ($1> 1) print}' | awk '{print $2}' | awk '$1 ~ /^UNRESOLVED/' > ~{input_vcf}_complex_events.bed

        #extract calls for the interesting events
        zcat ~{input_vcf}.bed.gz | grep -f ~{input_vcf}_complex_events > ~{input_vcf}_complex_events.bed

        cat ~{input_vcf}_complex_events.bed | awk -v OFS="\t"  '{print $1,$2,$3,$4,$6,$5}' > ~{input_vcf}_complex_events_formatted.bed

        rm ~{input_vcf}_complex_events.bed
        rm ~{input_vcf}.bed.gz
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

task clusterCPX{
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
        File svtk_bedcluster = "~{prefix}_cpx.gz"
    }
    command <<<
        set -euo pipefail

        ## combine all input beds into one for per batch allele count
        cat ~{input_beds} >> unified.bed
        svtk bedcluster unified.bed complex_unified_cluster.bed -f 0.5
        bgzip -c complex_unified_cluster.bed > complex_unified_cluster.bed.gz
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
