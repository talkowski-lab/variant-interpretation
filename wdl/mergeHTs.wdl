version 1.0

import "wes-denovo-helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow mergeHTs {
    input {
        Array[String] ht_uris
        String merged_filename
        String hail_docker
        Float? input_size
    }
    if (!defined(input_size)) {
        call helpers.getHailMTSizes as getHailHTSizes {
            input:
                mt_uris=ht_uris,
                hail_docker=hail_docker
        }
    }
    call helpers.mergeHTs as mergeHTs {
        input:
            ht_uris=ht_uris,
            merged_filename=merged_filename,
            hail_docker=hail_docker,
            input_size=select_first([input_size, getHailHTSizes.mt_size])
    }

    output {
        File merged_tsv = mergeHTs.merged_tsv
    }
}

