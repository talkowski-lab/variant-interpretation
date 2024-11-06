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
        String output_vcf_name
        # String hail_docker
        String sv_base_mini_docker
        Float input_size
    }

    call reheaderVCF_bcftools {
        input:
        new_header=new_header,
        vcf_uri=vcf_uri,
        output_vcf_name=output_vcf_name,
        sv_base_mini_docker=sv_base_mini_docker,
        input_size=input_size
    }

    output {
        File reheadered_vcf = reheaderVCF_bcftools.reheadered_vcf
    }
}

task reheaderVCF_bcftools {
    input {
        File new_header
        String vcf_uri
        String output_vcf_name
        String sv_base_mini_docker
        Float input_size
        RuntimeAttr? runtime_attr_override
    }
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
        mkfifo /tmp/token_fifo
        ( while true ; do curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token > /tmp/token_fifo ; done ) &
        HTS_AUTH_LOCATION=/tmp/token_fifo bcftools reheader -h ~{new_header} -o ~{output_vcf_name} ~{vcf_uri}
    >>>

    output {
        File reheadered_vcf = output_vcf_name
    }
}