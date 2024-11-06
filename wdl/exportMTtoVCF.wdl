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

workflow exportMTtoVCF {
    input {
        String mt_uri
        String hail_docker
    }

    call helpers.getHailMTSize as getHailMTSize {
        input:
            mt_uri=mt_uri,
            hail_docker=hail_docker
    }

    call exportMT {
        input:
            input_size=getHailMTSize.mt_size,
            mt_uri=mt_uri,
            hail_docker=hail_docker
    }

    output {
        File vcf_file = exportMT.vcf_file
    }
}

task exportMT {
    input {
        Float input_size
        String mt_uri
        String hail_docker
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
    cat <<EOF > export_mt.py
    import sys
    import hail as hl
    import numpy as np 
    import os

    mt_uri = sys.argv[1]
    cores = sys.argv[2]
    mem = int(np.floor(float(sys.argv[3])))

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")
        
    mt = hl.read_matrix_table(mt_uri)
    hl.export_vcf(mt, os.path.basename(mt_uri).split('.mt')[0]+'.vcf.bgz')
    EOF    

    python3 export_mt.py ~{mt_uri} ~{cpu_cores} ~{memory}
    >>>

    output {
        File vcf_file = basename(mt_uri, '.mt') + '.vcf.bgz'
    }
}