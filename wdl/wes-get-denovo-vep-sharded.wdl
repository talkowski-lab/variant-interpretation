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
        Array[String] de_novo_hts
        String merged_filename
        String hail_docker
    }

    Array[Pair[String, String]] mt_ht_pairs = zip(filtered_mts, de_novo_hts)

    scatter (mt_ht_pair in mt_ht_pairs) {
        String filtered_mt = mt_ht_pair.left
        String de_novo_ht = mt_ht_pair.right

        call getDenovoVEPandLOEUF.getDenovoVEP as getDenovoVEP {
            input:
            filtered_mt=filtered_mt,
            de_novo_ht=de_novo_ht,
            hail_docker=hail_docker
        }
    }

    call helpers.mergeResults as mergeResults {
        input:
        tsvs=getDenovoVEP.de_novo_vep,
        hail_docker=hail_docker,
        merged_filename=merged_filename
    }

    call getDenovoVEPandLOEUF.annotateLOEUF as annotateLOEUF {
        input:
        loeuf_file=loeuf_file,
        de_novo_vep=mergeResults.merged_tsv,
        hail_docker=hail_docker
    }

    output {
        File de_novo_vep = annotateLOEUF.de_novo_vep_loeuf
    }
}