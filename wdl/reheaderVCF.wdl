version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow reheaderVCF {
    input {
        File new_header
        File vcf_uri
        String output_vcf_uri
        String sv_base_mini_docker
    }

    call reheaderVCF_bcftools {
        input:
        new_header=new_header,
        vcf_uri=vcf_uri,
        output_vcf_uri=output_vcf_uri,
        sv_base_mini_docker=sv_base_mini_docker
    }
}

task reheaderVCF_bcftools {
    input {
        File new_header
        File vcf_uri
        String output_vcf_uri
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_uri, 'GB')
    Float base_disk_gb = 10.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size),
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -eou pipefail
        bcftools reheader -h ~{new_header} -o ~{basename(output_vcf_uri)} ~{vcf_uri}
        gsutil -m cp ~{basename(output_vcf_uri)} ~{output_vcf_uri}
    >>>
}