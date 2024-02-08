version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow compressHailMT {
    input {
        String mt_uri
        String hail_docker
        Float mt_size
    }

    call compressMT {
        input:
            mt_uri=mt_uri,
            hail_docker=hail_docker,
            mt_size=mt_size
    }

    output {
        File compressed_mt = compressMT.compressed_mt
    }
}

task compressMT {
    input {
        String mt_uri
        String hail_docker
        Float mt_size
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = mt_size
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
        gsutil -m cp -r ~{mt_uri} .
        tar -zcvf ~{basename(mt_uri)}.gz ~{basename(mt_uri)}
    >>>

    output {
        File compressed_mt = "~{basename(mt_uri)}.gz"
    }
}

