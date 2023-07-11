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
        File reference
        File reference_index
        String buffer
        String buffer_large
        String igv_docker
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_igv
    }

    call{

    }

    output{
        File tar_gz_pe = select_first(runIGV_evidence.pe_plots)
    }
}

task runIGV_evidence{
        input{
            File varfile
            String family
            File ped_file
            Array[String] samples
            File reference
            File reference_index
            String buffer
            String buffer_large
            String igv_docker
            String variant_interpretation_docker
            RuntimeAttr? runtime_attr_igv
        }

    Float input_size = size(select_all([varfile, ped_file]), "GB")
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
        File pe_plots="~{family}_pe_igv_plots.tar.gz"
        }
    }