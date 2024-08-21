version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow mergeSpliceAI {
    input {
        String snv_ht
        String indel_ht
        String output_ht_path
        String hail_docker

        Float input_size
    }

    call mergeSNVIndelSpliceAI {
        input:
        snv_ht=snv_ht,
        indel_ht=indel_ht,
        output_ht_path=output_ht_path,
        hail_docker=hail_docker,
        input_size=input_size
    }

    output {
        String merged_ht = output_ht_path
    }
}

task mergeSNVIndelSpliceAI {
    input {
        String snv_ht
        String indel_ht
        String output_ht_path
        String hail_docker
        Float input_size
        RuntimeAttr? runtime_attr_override
    }
    Float base_disk_gb = 10.0
    Float input_disk_scale = 10.0
    RuntimeAttr runtime_default = object {
        mem_gb: 8,
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
    set -eou pipefail
    cat <<EOF > merge.py
    from pyspark.sql import SparkSession
    import hail as hl
    import numpy as np
    import sys
    import ast
    import os

    snv_uri = sys.argv[1]
    indel_uri = sys.argv[2]
    output_ht_path = sys.argv[3]
    cores = sys.argv[4]  # string
    mem = int(np.floor(float(sys.argv[5])))

    # hl.init()
    hl.init(min_block_size=128, spark_conf={
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{int(np.floor(mem*0.8))}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    snv_ht = hl.read_table(spliceAI_snv_uri)
    indel_ht = hl.read_table(spliceAI_indel_uri)
    spliceAI_ht = snv_ht.union(indel_ht)

    spliceAI_ht.write(output_ht_path)
    EOF
    python3 merge.py ~{snv_ht} ~{indel_ht} ~{output_ht_path} ~{cpu_cores} ~{memory}
    >>>
}