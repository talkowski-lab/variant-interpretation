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
        File repetitive_regions_bed
        File rep_regions
        String var_type  # Indel or SNV
        String bagging_pu_source_script
        String bagging_pu_rf_script
        String tsv_to_bed_script
        String cohort_prefix
        String sv_base_mini_docker
        String hail_docker
        Array[String] metrics
        Array[Array[String]] features
        String known_vars_uri='false'
        Float vqslod_cutoff=-10
        Float prop_dn=1
        Int n_estimators_rf=100
        Int n_bag=10
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

    scatter (numeric in features) {
        scatter (metric in metrics) {
            call runBaggingPU_RF {
                input:
                    vcf_metrics_tsv_final=vcf_metrics_tsv_final,
                    ultra_rare_variants_tsv=ultra_rare_variants_tsv,
                    rep_regions=rep_regions,
                    ultra_rare_rep_regions=flagRepetitiveRegions.output_bed,
                    var_type=var_type,
                    bagging_pu_source_script=bagging_pu_source_script,
                    bagging_pu_rf_script=bagging_pu_rf_script,
                    cohort_prefix=cohort_prefix,
                    hail_docker=hail_docker,
                    metric=metric,
                    numeric=numeric,
                    vqslod_cutoff=vqslod_cutoff,
                    known_vars_uri=known_vars_uri,
                    prop_dn=prop_dn,
                    n_estimators_rf=n_estimators_rf,
                    n_bag=n_bag,
                    runtime_attr_override=runtime_attr_bagging_pu
            }
        }
    }
    
    output {
        Array[File] bagging_pu_results = flatten(runBaggingPU_RF.bagging_pu_results)
        Array[File] bagging_pu_importances = flatten(runBaggingPU_RF.bagging_pu_importances)
        Array[File] bagging_pu_oob_scores = flatten(runBaggingPU_RF.bagging_pu_oob_scores)
        # Array[File] bagging_pu_best_params = flatten(runBaggingPU_RF.bagging_pu_best_params)
        # Array[File] bagging_pu_figures = flatten(runBaggingPU_RF.bagging_pu_figures)
    }
}

task runBaggingPU_RF {
    input {
        File vcf_metrics_tsv_final
        File ultra_rare_variants_tsv
        File ultra_rare_rep_regions
        File rep_regions
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
        Int n_estimators_rf
        Int n_bag
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
        python3 run_bagging_pu.py ~{vcf_metrics_tsv_final} ~{ultra_rare_variants_tsv} \
        ~{cohort_prefix} ~{var_type} ~{sep=',' numeric} ~{vqslod_cutoff} \
        ~{prop_dn} ~{ultra_rare_rep_regions} ~{rep_regions} ~{known_vars_uri} ~{metric} \
        ~{n_estimators_rf} ~{n_bag} > stdout
    >>>

    output {
        File bagging_pu_results = "~{cohort_prefix}_~{var_type}_~{metric}_RF_results.tsv"
        File bagging_pu_importances = "~{cohort_prefix}_~{var_type}_~{metric}_RF_feature_importances.tsv"
        File bagging_pu_oob_scores = "~{cohort_prefix}_~{var_type}_~{metric}_RF_oob_scores.tsv"
        # File bagging_pu_best_params = "~{cohort_prefix}_~{var_type}_~{metric}_RF_best_params.tsv"
        Array[File] bagging_pu_figures = glob('*.png')
    }
}

# task compareMetrics {
#     input {
#         Array[File] bagging_pu_results
#         Array[String] metrics
#         String known_vars_uri
#         String hail_docker
#         RuntimeAttr? runtime_attr_override
#     }
    
#     Float input_size = size([bagging_pu_results], "GB")
#     Float base_disk_gb = 10.0
#     Float input_disk_scale = 5.0

#     RuntimeAttr runtime_default = object {
#         mem_gb: 4,
#         disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
#         cpu_cores: 1,
#         preemptible_tries: 3,
#         max_retries: 1,
#         boot_disk_gb: 10
#     }

#     RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

#     runtime {
#         memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
#         disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
#         cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
#         preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
#         maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
#         docker: hail_docker
#         bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
#     }

#     command <<<
#     cat <<EOF > compare_metrics.py
#     import pandas as pd
#     import seaborn as sns
#     import matplotlib.pyplot as plt
#     import sys

#     results_uri = sys.argv[1].split(',')
#     metrics = sys.argv[2].split(',')
#     cohort_prefix = sys.argv[3]
#     known_vars_uri = sys.argv[4]

#     known_vars_exist = (known_vars_uri!='false')

#     merged_results = pd.DataFrame()
#     for i, uri in enumerate(results_uri):
#         df = pd.read_csv(uri, sep='\t')
#         df['oob_score'] = metrics[i]
#         merged_results = pd.concat([merged_results, df])

#     fig, ax = plt.subplots(figsize=(8, 4));
#     if known_vars_exist:
#         sns.violinplot(merged_results,
#             x='oob_score', y='predict_proba_bag', hue='is_known', dodge=True, ax=ax);
#         ax.set_title(f"{cohort_prefix} Indels, Random Forest");
#     else:
#         sns.violinplot(merged_results,
#             x='oob_score', y='predict_proba_bag', hue='is_known', dodge=True, ax=ax);
#         ax.set_title(f"{cohort_prefix} Indels, Random Forest");
#     EOF
    
#     # python3 compare_metrics.py ~{}
#     >>>

#     output {
#         String temp = 'TODO'
#     }
# }