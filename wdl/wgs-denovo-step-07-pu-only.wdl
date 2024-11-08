version 1.0

# assuming filterUltraRareInheritedVariantsHail and filterUltraRareParentsVariantsHail already run
import  "flagRepetitiveRegions.wdl" as flagRepetitiveRegions
import "wgs-denovo-bagging-pu-rf-len.wdl" as BaggingPU_RF
import "filterUltraRareInheritedVariantsHail.wdl" as filterUltraRareInheritedVariantsHail
import "filterUltraRareParentsVariantsHail.wdl" as filterUltraRareParentsVariantsHail 

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow step7 {
    input {
        Array[File] annot_vcf_files
        File downsampled_ultra_rare_inherited
        File downsampled_ultra_rare_parents
        File lcr_uri
        File ped_sex_qc
        File meta_uri
        File trio_uri
        File vcf_metrics_tsv_final
        File hg38_reference
        File hg38_reference_dict
        File hg38_reference_fai
        String jvarkit_docker
        String hail_docker
        String sv_base_mini_docker
        String cohort_prefix

        RuntimeAttr? runtime_attr_filter_vcf
        RuntimeAttr? runtime_attr_merge_results
        RuntimeAttr? runtime_attr_downsample
    
        File repetitive_regions_bed
        String var_type  # Indel or SNV
        String bagging_pu_source_script
        String bagging_pu_rf_len_script
        String tsv_to_bed_script
        String cohort_prefix
        String metric='fp_fn_ratio'
        Array[String] sample_features=["GQ_parent", "AB_sample", "DPC_sample", "DPC_parent", "PL_sample_0.0", "PL_sample_1.1"]
        Array[String] variant_features=["MQ", "FS", "BaseQRankSum", "SOR", "LEN", "ReadPosRankSum", "DP", "QD", "VQSLOD"]
        Float vqslod_cutoff=-10
        Int n_estimators_rf=100
        Int n_bag=10
        Int n_jobs=-1
        Boolean filter_pass_before=false
        RuntimeAttr? runtime_attr_bagging_pu
    }

    call BaggingPU_RF.BaggingPU_RF as BaggingPU_RF {
        input:
        vcf_metrics_tsv_final=vcf_metrics_tsv_final,
        ultra_rare_inherited_tsv=downsampled_ultra_rare_inherited,
        ultra_rare_parents_tsv=downsampled_ultra_rare_parents,
        repetitive_regions_bed=repetitive_regions_bed,
        var_type=var_type,
        bagging_pu_source_script=bagging_pu_source_script,
        bagging_pu_rf_len_script=bagging_pu_rf_len_script,
        tsv_to_bed_script=tsv_to_bed_script,
        cohort_prefix=cohort_prefix,
        sv_base_mini_docker=sv_base_mini_docker,
        hail_docker=hail_docker,
        metric=metric,
        sample_features=sample_features,
        variant_features=variant_features,
        vqslod_cutoff=vqslod_cutoff,
        n_estimators_rf=n_estimators_rf,
        n_bag=n_bag,
        n_jobs=n_jobs,
        filter_pass_before=filter_pass_before,
        runtime_attr_bagging_pu=runtime_attr_bagging_pu
    }
    
    output {
        File vcf_metrics_tsv_final_pu = BaggingPU_RF.vcf_metrics_tsv_final_pu
    }
}