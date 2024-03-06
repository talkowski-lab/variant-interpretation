version 1.0

import "wes-denovo-helpers.wdl" as helpers
import "exportMTtoVCF.wdl" as exportMTtoVCF

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow exportMTtoVCF {
    input {
        Array[String] mt_uris
        String hail_docker
    }

    scatter (mt_uri in mt_uris) {
        call exportMTtoVCF.exportMTtoVCF as exportMT {
            input:
            mt_uri=mt_uri,
            hail_docker=hail_docker
        }
    }

    output {
        Array[File] vcf_files = exportMT.vcf_file
    }
}
