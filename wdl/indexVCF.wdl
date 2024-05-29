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
        Boolean localize_vcf
        Boolean use_tabix
        Boolean tbi
        RuntimeAttr? runtime_attr_index
    }

    String tbi_ext = if tbi then ".tbi" else ".csi"

    if (localize_vcf) {
        if (use_tabix) {
            call indexVCF_tabix {
                input:
                    vcf_uri=vcf_uri,
                    tbi=tbi,
                    sv_base_mini_docker=sv_base_mini_docker,
                    runtime_attr_override=runtime_attr_index
            }
        }

        if (!use_tabix) {
            call indexVCF_bcftools {
                input:
                    vcf_uri=vcf_uri,
                    tbi=tbi,
                    sv_base_mini_docker=sv_base_mini_docker,
                    runtime_attr_override=runtime_attr_index
            }
        }

        File vcf_idx_ = select_first([indexVCF_tabix.vcf_idx, indexVCF_bcftools.vcf_idx])
    }

    if (!localize_vcf) {
        if (use_tabix) {
            call indexVCF_tabix_remote {
                input:
                    vcf_uri=vcf_uri,
                    tbi=tbi,
                    sv_base_mini_docker=sv_base_mini_docker,
                    runtime_attr_override=runtime_attr_index
            }
        }

        if (!use_tabix) {
            call indexVCF_bcftools_remote {
                input:
                    vcf_uri=vcf_uri,
                    tbi=tbi,
                    sv_base_mini_docker=sv_base_mini_docker,
                    runtime_attr_override=runtime_attr_index
            }
        }
        String vcf_filename = vcf_uri
        File vcf_idx_remote = vcf_filename + tbi_ext
    }

    output {
        File vcf_idx = select_first([vcf_idx_, vcf_idx_remote])
    }
}

task indexVCF_tabix {
    input {
        File vcf_uri
        String sv_base_mini_docker
        Boolean tbi
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_uri, 'GB')
    Float base_disk_gb = 10.0
    Float base_mem_gb = 4.0
    Float input_disk_scale = 5.0
    
    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb,
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
    
    String tbi_str = if tbi then "" else "-C"
    String tbi_ext = if tbi then ".tbi" else ".csi"

    command <<<
        tabix ~{tbi_str} --verbosity 3 ~{vcf_uri}
    >>>

    output {
        File vcf_idx = basename(vcf_uri) + tbi_ext
    }
}


task indexVCF_bcftools {
    input {
        File vcf_uri
        String sv_base_mini_docker
        Boolean tbi
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_uri, 'GB')
    Float base_disk_gb = 10.0
    Float base_mem_gb = 4.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb,
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String tbi_str = if tbi then "-t" else ""
    String tbi_ext = if tbi then ".tbi" else ".csi"

    command <<<
        bcftools index ~{tbi_str} ~{vcf_uri}
    >>>

    output {
        File vcf_idx = basename(vcf_uri) + tbi_ext
    }
}

task indexVCF_tabix_remote {
    input {
        String vcf_uri
        String sv_base_mini_docker
        Boolean tbi
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
    
    String tbi_str = if tbi then "" else "-C"
    command <<<
        mkfifo /tmp/token_fifo
        ( while true ; do curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token > /tmp/token_fifo ; done ) &
        HTS_AUTH_LOCATION=/tmp/token_fifo tabix ~{tbi_str} --verbosity 3 ~{vcf_uri}
    >>>
}


task indexVCF_bcftools_remote {
    input {
        String vcf_uri
        String sv_base_mini_docker
        Boolean tbi
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

    String tbi_str = if tbi then "-t" else ""

    command <<<
        mkfifo /tmp/token_fifo
        ( while true ; do curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token > /tmp/token_fifo ; done ) &
        HTS_AUTH_LOCATION=/tmp/token_fifo bcftools index ~{tbi_str} ~{vcf_uri}
    >>>
}
