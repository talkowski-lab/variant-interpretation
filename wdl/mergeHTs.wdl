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
        String cohort_prefix
        String bucket_id
        String hail_docker
    }

    call helpers.mergeHTs as mergeHTs {
        input:
            ht_uris=ht_uris,
            cohort_prefix=cohort_prefix,
            bucket_id=bucket_id,
            hail_docker=hail_docker
    }

    output {
        String merged_ht = mergeHTs.merged_ht
    }
}

