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

workflow mergeTDT {
    input {
        Array[String] tdt_mt
        String cohort_prefix
        String hail_docker
        RuntimeAttr? runtime_attr_merge_tdt
    }

    call helpers.getHailMTSizes as getSize {
        input:
        mt_uris=tdt_mt,
        hail_docker=hail_docker
    }

    call mergeResultsTDT {
        input:
        tdt_mt=tdt_mt,
        input_size=getSize.mt_size,
        cohort_prefix=cohort_prefix,
        hail_docker=hail_docker,
        runtime_attr_override=runtime_attr_merge_tdt
    }

    output {
        File tdt_tsv_merged = mergeResultsTDT.tdt_tsv_merged
    }
}

task mergeResultsTDT {
    input {
        Array[String] tdt_mt
        String cohort_prefix
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
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
    cat << EOF > merge_tdt.py
        import hail as hl
        import pandas as pd
        import numpy as np
        import os
        import sys

        mt_list = sys.argv[1].split(',')
        cohort_prefix = sys.argv[2]

        tdt_df = pd.DataFrame()
        for i, uri in enumerate(mt_list):
            print(i)
            tdt_mt = hl.read_matrix_table(uri)
            tdt_df = pd.concat([tdt_df, tdt_mt.rows().to_pandas()])    
        tdt_df.to_csv(f"{cohort_prefix}_trio_tdt_merged.tsv", sep='\t', index=False)    
    EOF
    python3 merge_tdt.py ~{sep=',' tdt_mt} ~{cohort_prefix}
    >>>

    output {
        File tdt_tsv_merged = "~{cohort_prefix}_trio_tdt_merged.tsv"
    }
}