version 1.0

# IMPORT
## modified from ResolveCTX.wdl

import "Structs.wdl"
import "TinyResolve.wdl"

# WORKFLOW DEFINITION
workflow ResolveCTX {
    input {
        String prefix # batchid
        String docker_path
        Array[String] samples
        File manta_vcf_tar # batch std manta tarball
        File cytoband
        # File cytoband_idx
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
        RuntimeAttr? runtime_attr_cluster
    }

    call TinyResolve.TinyResolveCPX as TinyResolveCPX {
        input:
            samples = samples,
            manta_vcf_tar = manta_vcf_tar,
            cytoband = cytoband,
            # cytoband_idx = cytoband_idx,
            discfile = discfile,
            mei_bed = mei_bed,
            samples_per_shard = samples_per_shard,
            sv_pipeline_docker = sv_pipeline_docker,
            linux_docker = linux_docker,
            runtime_attr_resolve = runtime_attr_resolve,
            runtime_attr_untar = runtime_attr_untar
    }

    scatter (file in TinyResolveCPX.cpx_manta_unresolved_vcf) {
        String file_idx = basename(file, ".tbi")

        call extract_complex {
            input:
                input_vcf = file,
                docker = docker_path,
                runtime_attr_override = runtime_attr_extract
        }
    }

    #Array[File] all_cpx_formatted = extract_complex.cpx_formatted

    call clusterCPX {
        input:
            input_beds = extract_complex.cpx_formatted,
            input_dictionaries = extract_complex.cpx_dictionary,
            docker = docker_path,
            prefix = prefix,
            runtime_attr_override = runtime_attr_cluster
    }

    output {
        File cluster_bed = clusterCPX.svtk_bedcluster
        File cluster_dictionary = clusterCPX.svtk_dictionary
        File cluster_bed_filtered = clusterCPX.filtered_bedcluster
        Array[File] all_cpx_dictionary = extract_complex.cpx_dictionary
        Array[File] all_cpx_formatted = extract_complex.cpx_formatted
    }
}

# TASK DEFINITIONS

task extract_complex {
    input {
        File input_vcf
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(input_vcf, ".vcf.gz")

    RuntimeAttr default_attr = {
        "cpu_cores": 1,
        "mem_gb": 12,
        "disk_gb": 4,
        "boot_disk_gb": 8,
        "preemptible_tries": 3,
        "max_retries": 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File cpx_dictionary = "~{prefix}_complex_events.bed"
        File cpx_formatted = "~{prefix}_complex_events_formatted.bed"
    }

    command <<<
        set -euo pipefail

        # Convert to bed file
        svtk vcf2bed -i ALL --include-filters ~{input_vcf} ~{prefix}.bed
        bgzip -c ~{prefix}.bed > ~{prefix}.bed.gz

        # Extract multiple events
        zcat ~{prefix}.bed.gz | awk '{print $18}' | sort | uniq -c | \
        awk '{if ($1 > 1) print}' | awk '{print $2}' | awk '$1 ~ /^UNRESOLVED/' > ~{prefix}_complex_events

        # Extract calls for the interesting events
        zcat ~{prefix}.bed.gz | grep -f ~{prefix}_complex_events > ~{prefix}_complex_events.bed

        ## for svtk formatting, create +/-50bp padding + generate unique name that includes sample + event ID + make each name unique
        ## breakpoint 1
        cat ~{prefix}_complex_events.bed | awk -v OFS="\t" '{print $1"_"$7,$2-50,$3+50,$18"_"$4"_bp1",$6,$5}' > ~{prefix}_complex_events_formatted.bed1

        ## breakpoint2
        cat ~{prefix}_complex_events.bed | awk -v OFS="\t" '{print $8"_"$7,$9-50,$9+$11+50,$18"_"$4"_bp2",$6,$7}' > ~{prefix}_complex_events_formatted.bed2

        ## stitch
        cat ~{prefix}_complex_events_formatted.bed1 ~{prefix}_complex_events_formatted.bed2 > ~{prefix}_complex_events_formatted.bed
        
    >>>

    runtime {
        cpu: runtime_attr["cpu_cores"]
        memory: runtime_attr["mem_gb"] + " GiB"
        disks: "local-disk " + runtime_attr["disk_gb"] + " HDD"
        bootDiskSizeGb: runtime_attr["boot_disk_gb"]
        preemptible: runtime_attr["preemptible_tries"]
        maxRetries: runtime_attr["max_retries"]
        docker: docker
    }
}

task clusterCPX {
    input {
        Array[File] input_beds
        Array[File] input_dictionaries
        String docker
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = {
        "cpu_cores": 1,
        "mem_gb": 12,
        "disk_gb": 50,
        "boot_disk_gb": 8,
        "preemptible_tries": 3,
        "max_retries": 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File svtk_bedcluster = "~{prefix}_cpx.bed.gz"
        File svtk_dictionary = "~{prefix}_cpx.dictionary.bed.gz"
        File filtered_bedcluster = "~{prefix}_cpx_AC0.01.bed.gz"
    }

    command <<<
        set -euo pipefail

        ## combine all input beds into one for per batch allele count
        cat ~{sep=" " input_beds} > unified.bed
        svtk bedcluster unified.bed complex_unified_cluster.bed
        bgzip -c complex_unified_cluster.bed > ~{prefix}_cpx.bed.gz

        ## generate per batch dictionary for post filtering
        cat ~{sep=" " input_dictionaries} > dictionary.bed
        bgzip -c dictionary.bed > ~{prefix}_cpx.dictionary.bed.gz
        Rscript /opt/cpx_AC.R
        bgzip -c filtered.bed > ~{prefix}_cpx_AC0.01.bed.gz

    >>>

    runtime {
        cpu: runtime_attr["cpu_cores"]
        memory: runtime_attr["mem_gb"] + " GiB"
        disks: "local-disk " + runtime_attr["disk_gb"] + " HDD"
        bootDiskSizeGb: runtime_attr["boot_disk_gb"]
        preemptible: runtime_attr["preemptible_tries"]
        maxRetries: runtime_attr["max_retries"]
        docker: docker
    }
}
