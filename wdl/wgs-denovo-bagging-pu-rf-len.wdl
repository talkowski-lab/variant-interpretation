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
        File ultra_rare_inherited_tsv
        File ultra_rare_parents_tsv
        File repetitive_regions_bed
        String var_type  # Indel or SNV
        String bagging_pu_source_script
        String bagging_pu_rf_len_script
        String tsv_to_bed_script
        String cohort_prefix
        String sv_base_mini_docker
        String hail_docker
        String metric
        Array[String] sample_features
        Array[String] variant_features
        Float vqslod_cutoff=-10
        Float prop_dn=1
        Int n_estimators_rf=100
        Int n_bag=10
        Boolean filter_pass_before=false
        RuntimeAttr? runtime_attr_bagging_pu
    }

    call flagRepetitiveRegions.flagRepetitiveRegions as flagRepetitiveRegions {
        input:
        tsv=vcf_metrics_tsv_final,
        repetitive_regions_bed=repetitive_regions_bed,
        tsv_to_bed_script=tsv_to_bed_script,
        sv_base_mini_docker=sv_base_mini_docker,
        hail_docker=hail_docker
    }

    call runBaggingPU_RF {
        input:
            vcf_metrics_tsv_final=vcf_metrics_tsv_final,
            ultra_rare_inherited_tsv=ultra_rare_inherited_tsv,
            ultra_rare_parents_tsv=ultra_rare_parents_tsv,
            rep_regions=flagRepetitiveRegions.output_bed,
            var_type=var_type,
            bagging_pu_source_script=bagging_pu_source_script,
            bagging_pu_rf_len_script=bagging_pu_rf_len_script,
            cohort_prefix=cohort_prefix,
            hail_docker=hail_docker,
            metric=metric,
            variant_features=variant_features,
            sample_features=sample_features,
            vqslod_cutoff=vqslod_cutoff,
            prop_dn=prop_dn,
            n_estimators_rf=n_estimators_rf,
            n_bag=n_bag,
            filter_pass_before=filter_pass_before,
            runtime_attr_override=runtime_attr_bagging_pu
    }
    
    output {
        File vcf_metrics_tsv_final_pu = runBaggingPU_RF.vcf_metrics_tsv_final_pu
    }
}

task runBaggingPU_RF {
    input {
        File vcf_metrics_tsv_final
        File ultra_rare_inherited_tsv
        File ultra_rare_parents_tsv
        File rep_regions
        Array[String] variant_features
        Array[String] sample_features
        String var_type
        String bagging_pu_source_script
        String bagging_pu_rf_len_script
        String cohort_prefix
        String hail_docker
        String metric
        Float vqslod_cutoff
        Float prop_dn
        Int n_estimators_rf
        Int n_bag
        Boolean filter_pass_before
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([vcf_metrics_tsv_final, ultra_rare_inherited_tsv], "GB")
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
        curl ~{bagging_pu_rf_len_script} > run_bagging_pu.py
        curl ~{bagging_pu_source_script} > baggingPU.py
        python3 run_bagging_pu.py ~{vcf_metrics_tsv_final} ~{ultra_rare_inherited_tsv} \
        ~{cohort_prefix} ~{var_type} ~{sep=',' variant_features} ~{sep=',' sample_features} ~{vqslod_cutoff} \
        ~{prop_dn} ~{rep_regions} ~{metric} \
        ~{n_estimators_rf} ~{n_bag} ~{filter_pass_before} > stdout
    >>>

    output {
        File vcf_metrics_tsv_final_pu = "~{basename(vcf_metrics_tsv_final, '.tsv.gz')}_pu_~{var_type}.tsv"
    }
}