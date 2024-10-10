version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow prioritizeCSQ {
    input {
        File vcf_metrics_tsv
        File vep_vcf_file
        String prioritize_csq_script
        String hail_docker
        String sample_column
        String genome_build
    }
    
    call annotateMostSevereCSQ {
        input:
        vcf_metrics_tsv=vcf_metrics_tsv,
        vep_vcf_file=vep_vcf_file,
        prioritize_csq_script=prioritize_csq_script,
        hail_docker=hail_docker,
        sample_column=sample_column,
        genome_build=genome_build
    }

    output {
        File vcf_metrics_tsv_prior_csq = annotateMostSevereCSQ.vcf_metrics_tsv_prior_csq
    }
}

task annotateMostSevereCSQ {
    input {
        File vcf_metrics_tsv
        File vep_vcf_file
        String prioritize_csq_script
        String hail_docker
        String sample_column
        String genome_build
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_metrics_tsv, "GB")
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
    
    String file_ext = if sub(basename(vcf_metrics_tsv), '\\.gz', '')==basename(vcf_metrics_tsv) then '.tsv' else '.tsv.gz'
    command <<<
        set -eou pipefail
        curl ~{prioritize_csq_script} > prioritize_csq.py
        python3 prioritize_csq.py ~{vcf_metrics_tsv} ~{cpu_cores} ~{memory} \
        ~{sample_column} ~{vep_vcf_file} ~{genome_build}
    >>>

    output {
        File vcf_metrics_tsv_prior_csq = basename(vcf_metrics_tsv, file_ext) + '_prioritized_csq.tsv.gz'
    }
}