version 1.0

import  "flagRepetitiveRegions.wdl" as flagRepetitiveRegions

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow BaggingPU_RF {
    input {
        File vcf_metrics_tsv_final
        File ultra_rare_variants_tsv
        File ultra_rare_polyx_vcf
        File repetitive_regions_bed
        String var_type
        String bagging_pu_source_script
        String bagging_pu_rf_script
        String tsv_to_bed_script
        String cohort_prefix
        String sv_base_mini_docker
        String hail_docker
        String metric
        Array[String]? numeric
        String known_vars_uri='false'
        Float vqslod_cutoff=-10
        Float prop_dn=1
        RuntimeAttr? runtime_attr_bagging_pu
    }

    call flagRepetitiveRegions.flagRepetitiveRegions as flagRepetitiveRegions {
        input:
        tsv=ultra_rare_variants_tsv,
        repetitive_regions_bed=repetitive_regions_bed,
        tsv_to_bed_script=tsv_to_bed_script,
        sv_base_mini_docker=sv_base_mini_docker,
        hail_docker=hail_docker
    }

    if (!defined(numeric)) {
        Array[String] numeric_default = ['false']
    }
    Array[String] numeric_ = select_first([numeric, numeric_default])

    call runBaggingPU_RF {
        input:
            vcf_metrics_tsv_final=vcf_metrics_tsv_final,
            ultra_rare_variants_tsv=ultra_rare_variants_tsv,
            ultra_rare_polyx_vcf=ultra_rare_polyx_vcf,
            repetitive_regions_bed=repetitive_regions_bed,
            ultra_rare_repetitive_regions_bed=flagRepetitiveRegions.output_bed,
            var_type=var_type,
            bagging_pu_source_script=bagging_pu_source_script,
            bagging_pu_rf_script=bagging_pu_rf_script,
            cohort_prefix=cohort_prefix,
            hail_docker=hail_docker,
            metric=metric,
            numeric=numeric_,
            vqslod_cutoff=vqslod_cutoff,
            known_vars_uri=known_vars_uri,
            prop_dn=prop_dn,
            runtime_attr_override=runtime_attr_bagging_pu
    }

    output {
        File bagging_pu_results = runBaggingPU_RF.bagging_pu_results
        File bagging_pu_importances = runBaggingPU_RF.bagging_pu_importances
        File bagging_pu_oob_scores = runBaggingPU_RF.bagging_pu_oob_scores
        File bagging_pu_best_params = runBaggingPU_RF.bagging_pu_best_params
        Array[File] bagging_pu_figures = runBaggingPU_RF.bagging_pu_figures
    }
}

task runBaggingPU_RF {
    input {
        File vcf_metrics_tsv_final
        File ultra_rare_variants_tsv
        File ultra_rare_polyx_vcf
        File ultra_rare_repetitive_regions_bed
        File repetitive_regions_bed
        Array[String] numeric
        String var_type
        String bagging_pu_source_script
        String bagging_pu_rf_script
        String cohort_prefix
        String hail_docker
        String metric
        String known_vars_uri
        Float vqslod_cutoff
        Float prop_dn
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([vcf_metrics_tsv_final, ultra_rare_variants_tsv], "GB")
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
        curl ~{bagging_pu_rf_script} > run_bagging_pu.py
        curl ~{bagging_pu_source_script} > baggingPU.py
        python3 run_bagging_pu.py ~{vcf_metrics_tsv_final} ~{ultra_rare_variants_tsv} ~{ultra_rare_polyx_vcf} \
        ~{cohort_prefix} ~{var_type} ~{sep=',' numeric} ~{vqslod_cutoff} \
        ~{prop_dn} ~{ultra_rare_repetitive_regions_bed} ~{repetitive_regions_bed} ~{known_vars_uri} ~{metric}
    >>>

    output {
        File bagging_pu_results = "~{cohort_prefix}_~{var_type}_RF_results.tsv"
        File bagging_pu_importances = "~{cohort_prefix}_~{var_type}_RF_feature_importances.tsv"
        File bagging_pu_oob_scores = "~{cohort_prefix}_~{var_type}_RF_oob_scores.tsv"
        File bagging_pu_best_params = "~{cohort_prefix}_~{var_type}_RF_best_params.tsv"
        Array[File] bagging_pu_figures = glob('*.png')
    }
}

