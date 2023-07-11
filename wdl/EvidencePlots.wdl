#Author: Alba Sanchis-Juan <asanchis@broadinstitute.org>
#This script makes SV plots from PE SR and RD evidence


version 1.0

import "Structs.wdl"

workflow IGV_evidence {
    input{
        File varfile
        String family
        File ped_file
        Array[String] samples
        Array[File] disc_files
        Array[File] split_files
        File reference
        File reference_index
        String buffer
        String buffer_large
        String igv_docker
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_reformat_pesr
    }

    call reformatPESR{
        input:
            varfile=varfile,
            disc_files=disc_files,
            split_files=split_files,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override=runtime_attr_reformat_pesr

    }

    output{
        File tar_gz_pe = select_first(runIGV_evidence.pe_plots)
    }
}

task reformatPE{
    input{
        File varfile
        File disc_file
        File split_file
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([varfile, disc_files, split_files]), "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + input_size),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
            set -euo pipefail

            ##Reformat disc reads
            zcat ~{file} | \
                awk '$1 == $4 {print $1"\t"$2"\t"$5"\t.\t1\t"$3"\t255,0,0"}' | \
                sed -e 's/-\t255,0,0/-\t0,255,0/g' | \
                bgzip -c >  tmp.bed.gz

            tabix -p bed tmp.bed.gz
            bedtools intersect -a tmp.bed.gz -b ~{varfile} -f 0.5 -wa > pe_roi.junctions.bed
            >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: igv_docker
    }
    output{
        File pe_reformat="pe_roi.junctions.bed"
        }
}


task reformatSR{
    input{
        File varfile
        File disc_file
        File split_file
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([varfile, disc_files, split_files]), "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + input_size),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
            set -euo pipefail

            ##Reformat split reads
            zcat ~{file} | awk -F"\t" 'BEGIN { OFS="\t" }{print $1,$2-5,$2+5,$3,$4}' | \
                sed -e 's/left/-/g' | \
                sed -e 's/right/+/g' | \
                bgzip -c > tmp.bed.gz

            tabix -p bed tmp.bed.gz
            tabix -R ~{varfile} tmp.bed.gz | \
                bgzip -c > sr_roi.junctions.bed
            >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: igv_docker
    }
    output{
        File sr_evidence="sr_roi.junctions.bed"
        }
}