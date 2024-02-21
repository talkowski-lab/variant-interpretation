version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow runBaggingPU {
    input {
        File vcf_metrics_tsv
        File ultra_rare_variants_tsv
        Array[String] outlier_samples
        String var_type
        String bagging_pu_source_script
        String run_bagging_pu_script
        String cohort_prefix
        String hail_docker
        Array[String]? numeric
        RuntimeAttr? runtime_attr_bagging_pu
    }

    if (!defined(numeric)) {
        Array[String] numeric_default = ['false']
    }
    Array[String] numeric_ = select_first([numeric, numeric_default])

    call baggingPU {
        input:
            vcf_metrics_tsv=vcf_metrics_tsv,
            ultra_rare_variants_tsv=ultra_rare_variants_tsv,
            outlier_samples=outlier_samples,
            var_type=var_type,
            bagging_pu_source_script=bagging_pu_source_script,
            run_bagging_pu_script=run_bagging_pu_script,
            cohort_prefix=cohort_prefix,
            hail_docker=hail_docker,
            numeric=numeric_,
            runtime_attr_override=runtime_attr_bagging_pu
    }

    output {
        File bagging_pu_results = baggingPU.bagging_pu_results
        File bagging_pu_importances = baggingPU.bagging_pu_importances
    }
}

task baggingPU {
    input {
        File vcf_metrics_tsv
        File ultra_rare_variants_tsv
        Array[String] outlier_samples
        String var_type
        String bagging_pu_source_script
        String run_bagging_pu_script
        String cohort_prefix
        String hail_docker
        Array[String] numeric
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([vcf_metrics_tsv, ultra_rare_variants_tsv], "GB")
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

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
    
    command <<<
        curl ~{run_bagging_pu_script} > run_bagging_pu.py
        curl ~{bagging_pu_source_script} > baggingPU.py
        python3 run_bagging_pu.py ~{vcf_metrics_tsv} ~{ultra_rare_variants_tsv} ~{cohort_prefix} \
        ~{var_type} ~{sep=',' numeric} ~{sep=',' outlier_samples}
    >>>

    output {
        File bagging_pu_results = "~{cohort_prefix}_baggingPU_~{var_type}_results.tsv"
        File bagging_pu_importances = "~{cohort_prefix}_~{var_type}_feature_importances.tsv"
    }
}