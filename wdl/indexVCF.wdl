version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow indexVCF {
    input {
        String vcf_uri
        String sv_base_mini_docker
        Boolean use_tabix
        RuntimeAttr? runtime_attr_index
    }
    if (use_tabix) {
        call indexVCF_tabix {
            input:
                vcf_uri=vcf_uri,
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_index
        }
    }

    if (!use_tabix) {
    call indexVCF_bcftools {
        input:
            vcf_uri=vcf_uri,
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_index
    }
}

}

task indexVCF_tabix {
    input {
        String vcf_uri
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }
    
    Float base_disk_gb = 10.0
    Float base_mem_gb = 4.0
    
    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb,
        disk_gb: ceil(base_disk_gb),
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        mkfifo /tmp/token_fifo
        ( while true ; do curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token > /tmp/token_fifo ; done ) &
        HTS_AUTH_LOCATION=/tmp/token_fifo tabix --verbosity 3 ~{vcf_uri}
    >>>
}


task indexVCF_bcftools {
    input {
        String vcf_uri
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }
    
    Float base_disk_gb = 10.0
    Float base_mem_gb = 4.0
    
    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb,
        disk_gb: ceil(base_disk_gb),
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        mkfifo /tmp/token_fifo
        ( while true ; do curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token > /tmp/token_fifo ; done ) &
        HTS_AUTH_LOCATION=/tmp/token_fifo bcftools index -t ~{vcf_uri}
    >>>
}
