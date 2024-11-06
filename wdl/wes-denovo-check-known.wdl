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

workflow DenovoCheckKnown {
    input {
        String known_vars_uri
        String cohort_prefix
        Array[String] annot_mt
        Array[String] filtered_mt
        String hail_docker
        Float? input_size
    }

    scatter (annot_mt_uri in annot_mt) {
        if (!defined(input_size)) {
            call helpers.getHailMTSize as getAnnotMTSize {
                input:
                mt_uri=annot_mt_uri,
                hail_docker=hail_docker
            }
        }
        call overlapKnownMT as overlapAnnotMT {
            input:
            mt_uri=annot_mt_uri,
            known_vars_uri=known_vars_uri,
            hail_docker=hail_docker,
            input_size=select_first([input_size, getAnnotMTSize.mt_size])
        }
    }

    scatter (filtered_mt_uri in filtered_mt) {
        if (!defined(input_size)) {
            call helpers.getHailMTSize as getFilteredMTSize {
                input:
                mt_uri=filtered_mt_uri,
                hail_docker=hail_docker
            }
        }
        call overlapKnownMT as overlapFilteredMT {
            input:
            mt_uri=filtered_mt_uri,
            known_vars_uri=known_vars_uri,
            hail_docker=hail_docker,
            input_size=select_first([input_size, getFilteredMTSize.mt_size])
        }
    }

    call helpers.getHailMTSizes as getTotalAnnotMTSize {
        input:
        mt_uris=overlapAnnotMT.known_tsv,
        hail_docker=hail_docker
    }
    call helpers.getHailMTSizes as getTotalFilteredMTSize {
        input:
        mt_uris=overlapFilteredMT.known_tsv,
        hail_docker=hail_docker
    }

    call helpers.mergeResultsPython as mergeAnnotMT {
        input:
        tsvs=overlapAnnotMT.known_tsv,
        hail_docker=hail_docker,
        merged_filename=cohort_prefix + '_annot_mt_known_variants.tsv',
        input_size=getTotalAnnotMTSize.mt_size
    }
    call helpers.mergeResultsPython as mergeFilteredMT {
        input:
        tsvs=overlapFilteredMT.known_tsv,
        hail_docker=hail_docker,
        merged_filename=cohort_prefix + '_filtered_mt_known_variants.tsv',
        input_size=getTotalFilteredMTSize.mt_size
    }

    output {
        File annot_mt_known = mergeAnnotMT.merged_tsv
        File filtered_mt_known = mergeFilteredMT.merged_tsv
    }
}

task overlapKnownMT {
    input {
        String mt_uri
        String known_vars_uri
        String hail_docker
        Float input_size
        RuntimeAttr? runtime_attr_override
    }
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
    cat <<EOF > overlap_mt.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    mt_uri = sys.argv[1]
    known_vars_uri = sys.argv[2]
    cores = sys.argv[3]
    mem = int(np.floor(float(sys.argv[4])))

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    known_ht = hl.read_table(known_vars_uri)
    mt = hl.read_matrix_table(mt_uri)

    mt = mt.semi_join_rows(known_ht)
    df = mt.rows().to_pandas()
    filename = os.path.basename(mt_uri).split('.mt')[0] + '_known_variants.tsv'
    df.to_csv(filename, sep='\t', index=False)
    EOF

    python3 overlap_mt.py ~{mt_uri} ~{known_vars_uri} ~{cpu_cores} ~{memory} > stdout
    >>>

    output {
        File known_tsv = basename(mt_uri, '.mt') + '_known_variants.tsv'
    }
}