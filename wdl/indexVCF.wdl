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
        File vcf
        String vcf_dir
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_index
    }

    call indexVCF_ {
        input:
            vcf=vcf,
            vcf_dir=vcf_dir,
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_index
    }
}

task indexVCF_ {
    input {
        File vcf
        String vcf_dir
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }
    
    Float input_size = size(vcf, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 2.0
    Float input_mem_scale = 3.0
    Float input_disk_scale = 5.0
    
    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb + input_size * input_mem_scale,
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

    command {
        tabix ~{vcf}
        gsutil -m cp ~{vcf} ~{vcf_dir}/
        gsutil -m cp ~{vcf}.tbi ~{vcf_dir}/
    }
}
