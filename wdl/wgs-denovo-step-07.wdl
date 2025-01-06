version 1.0

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
        File lcr_uri
        File ped_sex_qc
        File meta_uri
        File trio_uri
        File vcf_metrics_tsv_final
        File hg38_reference
        File hg38_reference_dict
        File hg38_reference_fai
        File remove_regions_bed
        String prioritize_csq_script
        String filter_rare_inherited_python_script
        String filter_rare_parents_python_script
        String jvarkit_docker
        String hail_docker
        String sv_base_mini_docker
        String cohort_prefix
        Float AF_threshold=0.005
        Int AC_threshold=2
        Float csq_af_threshold=0.01
        Int gq_het_threshold=99
        Int gq_hom_ref_threshold=30
        Int qual_threshold=150
        Float sor_threshold_indel=3.0
        Float sor_threshold_snv=2.5
        Float readposranksum_threshold_indel=-1.7
        Float readposranksum_threshold_snv=-1.4
        Float qd_threshold_indel=4.0
        Float qd_threshold_snv=3.0
        Float mq_threshold=50
        Int shards_per_chunk=10

        # for downsampling
        Int chunk_size=100000
        Float snv_scale=1
        Float indel_scale=1
        # Boolean prioritize_gnomad=false

        RuntimeAttr? runtime_attr_filter_vcf
        RuntimeAttr? runtime_attr_merge_results
        RuntimeAttr? runtime_attr_downsample
    
        File repetitive_regions_bed
        # String var_type  # Indel or SNV, only Indel for now
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
        Boolean filter_pass_before=false
        RuntimeAttr? runtime_attr_bagging_pu
    }

    call filterUltraRareInheritedVariantsHail.filterUltraRareInheritedVariantsHail as filterUltraRareInheritedVariantsHail {
        input:
            annot_vcf_files=annot_vcf_files,
            lcr_uri=lcr_uri,
            ped_sex_qc=ped_sex_qc,
            meta_uri=meta_uri,
            trio_uri=trio_uri,
            vcf_metrics_tsv_final=vcf_metrics_tsv_final,
            hg38_reference=hg38_reference,
            hg38_reference_dict=hg38_reference_dict,
            hg38_reference_fai=hg38_reference_fai,
            remove_regions_bed=remove_regions_bed,
            prioritize_csq_script=prioritize_csq_script,
            filter_rare_inherited_python_script=filter_rare_inherited_python_script,
            jvarkit_docker=jvarkit_docker,
            hail_docker=hail_docker,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix,
            AF_threshold=AF_threshold,
            AC_threshold=AC_threshold,
            csq_af_threshold=csq_af_threshold,
            gq_het_threshold=gq_het_threshold,
            gq_hom_ref_threshold=gq_hom_ref_threshold,
            qual_threshold=qual_threshold,
            sor_threshold_indel=sor_threshold_indel,
            sor_threshold_snv=sor_threshold_snv,
            readposranksum_threshold_indel=readposranksum_threshold_indel,
            readposranksum_threshold_snv=readposranksum_threshold_snv,
            qd_threshold_indel=qd_threshold_indel,
            qd_threshold_snv=qd_threshold_snv,
            mq_threshold=mq_threshold,
            shards_per_chunk=shards_per_chunk,
            chunk_size=chunk_size,
            snv_scale=snv_scale,
            indel_scale=indel_scale,
            prioritize_gnomad=false,
            runtime_attr_filter_vcf=runtime_attr_filter_vcf,
            runtime_attr_merge_results=runtime_attr_merge_results,
            runtime_attr_downsample=runtime_attr_downsample
    }

    call filterUltraRareParentsVariantsHail.filterUltraRareParentsVariantsHail as filterUltraRareParentsVariantsHail {
        input:
            annot_vcf_files=annot_vcf_files,
            lcr_uri=lcr_uri,
            ped_sex_qc=ped_sex_qc,
            meta_uri=meta_uri,
            trio_uri=trio_uri,
            vcf_metrics_tsv_final=vcf_metrics_tsv_final,
            hg38_reference=hg38_reference,
            hg38_reference_dict=hg38_reference_dict,
            hg38_reference_fai=hg38_reference_fai,
            remove_regions_bed=remove_regions_bed,
            prioritize_csq_script=prioritize_csq_script,
            filter_rare_parents_python_script=filter_rare_parents_python_script,
            jvarkit_docker=jvarkit_docker,
            hail_docker=hail_docker,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix,
            AF_threshold=AF_threshold,
            AC_threshold=AC_threshold,
            csq_af_threshold=csq_af_threshold,
            gq_het_threshold=gq_het_threshold,
            gq_hom_ref_threshold=gq_hom_ref_threshold,
            qual_threshold=qual_threshold,
            sor_threshold_indel=sor_threshold_indel,
            sor_threshold_snv=sor_threshold_snv,
            readposranksum_threshold_indel=readposranksum_threshold_indel,
            readposranksum_threshold_snv=readposranksum_threshold_snv,
            qd_threshold_indel=qd_threshold_indel,
            qd_threshold_snv=qd_threshold_snv,
            mq_threshold=mq_threshold,
            shards_per_chunk=shards_per_chunk,
            chunk_size=chunk_size,
            snv_scale=snv_scale,
            indel_scale=indel_scale,
            prioritize_gnomad=true,
            runtime_attr_filter_vcf=runtime_attr_filter_vcf,
            runtime_attr_merge_results=runtime_attr_merge_results,
            runtime_attr_downsample=runtime_attr_downsample
    }

    call BaggingPU_RF.BaggingPU_RF as BaggingPU_RF {
        input:
        vcf_metrics_tsv_final=vcf_metrics_tsv_final,
        ultra_rare_inherited_tsv=filterUltraRareInheritedVariantsHail.downsampled_ultra_rare_inherited_Indel,
        ultra_rare_parents_tsv=filterUltraRareParentsVariantsHail.downsampled_ultra_rare_parents_Indel,
        repetitive_regions_bed=repetitive_regions_bed,
        var_type='Indel',  # Indel or SNV, only Indel for now
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
        filter_pass_before=filter_pass_before,
        runtime_attr_bagging_pu=runtime_attr_bagging_pu
    }
    
    output {
        File ultra_rare_inherited_tsv = filterUltraRareInheritedVariantsHail.ultra_rare_inherited_tsv
        File downsampled_ultra_rare_inherited_SNV = filterUltraRareInheritedVariantsHail.downsampled_ultra_rare_inherited_SNV
        File downsampled_ultra_rare_inherited_Indel = filterUltraRareInheritedVariantsHail.downsampled_ultra_rare_inherited_Indel

        File ultra_rare_parents_tsv = filterUltraRareParentsVariantsHail.ultra_rare_parents_tsv
        File downsampled_ultra_rare_parents_SNV = filterUltraRareParentsVariantsHail.downsampled_ultra_rare_parents_SNV
        File downsampled_ultra_rare_parents_Indel = filterUltraRareParentsVariantsHail.downsampled_ultra_rare_parents_Indel
        
        File vcf_metrics_tsv_final_pu = BaggingPU_RF.vcf_metrics_tsv_final_pu
    }
}