version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow filterRareVariantsHail {
    input {
        File vcf_file
        File lcr_uri
        File ped_uri
        File meta_uri
        File trio_uri
        File filter_rare_variants_python_script
        String vep_hail_docker
        String cohort_prefix
        RuntimeAttr? runtime_attr_filter_vcf
    }

    call filterRareVariants {
        input:
            vcf_file=vcf_file,
            lcr_uri=lcr_uri,
            ped_uri=ped_uri,
            meta_uri=meta_uri,
            trio_uri=trio_uri,
            filter_rare_variants_python_script=filter_rare_variants_python_script,
            vep_hail_docker=vep_hail_docker,
            cohort_prefix=cohort_prefix,
            runtime_attr_override=runtime_attr_filter_vcf
    }

    output {
        File ultra_rare_variants_tsv = filterRareVariants.ultra_rare_variants_tsv
    }
}

task filterRareVariants {
    input {
        File vcf_file
        File lcr_uri
        File ped_uri
        File meta_uri
        File trio_uri
        File filter_rare_variants_python_script
        String vep_hail_docker
        String cohort_prefix
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
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
        docker: vep_hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command {
        python 3.9 ~{filter_rare_variants_python_script} ~{lcr_uri} ~{ped_uri} ~{meta_uri} ~{trio_uri} ~{vcf_file} \
        ~{cohort_prefix} ~{cpu_cores} ~{memory}
    }

    output {
        File ultra_rare_variants_tsv = cohort_prefix + '_ultra_rare_variants.tsv'
    }
}