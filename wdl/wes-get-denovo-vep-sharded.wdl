version 1.0

import "wes-denovo-helpers.wdl" as helpers
import "wes-get-denovo-vep.wdl" as getDenovoVEPandLOEUF

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow getDenovoVEPandLOEUFSharded {
    input {
        File loeuf_file
        Array[String] filtered_mts
        Array[String] denovo_hts
        String merged_filename
        String hail_docker
    }

    Array[Pair[String, String]] mt_ht_pairs = zip(filtered_mts, denovo_hts)

    scatter (mt_ht_pair in mt_ht_pairs) {
        String filtered_mt = mt_ht_pair.left
        String denovo_ht = mt_ht_pair.right

        call getDenovoVEPandLOEUF.getDenovoVEPandLOEUF as getDenovoVEPandLOEUF {
            input:
            filtered_mt=filtered_mt,
            denovo_ht=denovo_ht,
            loeuf_file=loeuf_file,
            hail_docker=hail_docker
        }
    }

    call helpers.mergeResults as mergeResults {
        input:
        tsvs=getDenovoVEPandLOEUF.de_novo_vep,
        hail_docker=hail_docker,
        merged_filename=merged_filename
    }

    output {
        File de_novo_vep = mergeResults.merged_tsv
    }
}