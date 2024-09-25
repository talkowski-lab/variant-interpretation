version 1.0

import "wgs-denovo-step-01.wdl" as step1
import "wgs-denovo-step-02.wdl" as step2
import "wgs-denovo-step-03.wdl" as step3
import "wgs-denovo-step-04.wdl" as step4
import "wgs-denovo-step-05.wdl" as step5
import "wgs-denovo-step-06.wdl" as step6
import "wgs-denovo-step-07.wdl" as step7
import "annotateHPandVAF.wdl" as annotateHPandVAF
import "annotateMPCandLOEUF.wdl" as annotateMPCandLOEUF

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow wgs_denovo_full {
    input {
        String python_trio_sample_script
        String python_preprocess_script
        String subset_ped_script
        String uberSplit_v3_script
        String merge_vcf_to_tsv_fullQC_script
        String get_sample_pedigree_script
        String filter_final_tsv_script
        String annotate_mpc_loeuf_script
        String filter_rare_inherited_python_script
        String filter_rare_parents_python_script
        String bagging_pu_source_script
        String bagging_pu_rf_len_script
        String tsv_to_bed_script

        String mpc_dir
        File mpc_chr22_file
        File loeuf_file
        File lcr_uri
        File ped_sex_qc
        File relatedness_qc
        File hg38_reference
        File hg38_reference_fai
        File hg38_reference_dict
        File info_header
        File repetitive_regions_bed

        Array[File]? vep_vcf_files
        Array[File]? vep_annotated_final_vcf
        String cohort_prefix
        String sv_base_mini_docker
        String trio_denovo_docker
        String hail_docker
        String jvarkit_docker
        String sample_column
        Int batch_size

        Boolean exclude_gq_filters=false
        Boolean merge_split_vcf=false
        String genome_build='GRCh38'
        Int shards_per_chunk=10
        Int qual_threshold=150
        Float sor_threshold_indel=3.0
        Float sor_threshold_snv=2.5
        Float readposranksum_threshold_indel=-1.7
        Float readposranksum_threshold_snv=-1.4
        Float qd_threshold_indel=4.0
        Float qd_threshold_snv=3.0
        Float mq_threshold=50
        Float minDQ
        Float AF_threshold=0.005
        Int AC_threshold=2
        Float csq_af_threshold=0.01

        # for ultra-rare
        Int gq_het_threshold=99
        Int gq_hom_ref_threshold=30
    
        # String var_type  # Indel or SNV, only Indel for now
        String metric
        Array[String] sample_features
        Array[String] variant_features
        Float vqslod_cutoff=-10
        Int n_estimators_rf=100
        Int n_bag=10
        Boolean filter_pass_before=false
        RuntimeAttr? runtime_attr_bagging_pu
        RuntimeAttr? runtime_attr_merge_chunk
        RuntimeAttr? runtime_attr_filter_vcf
        RuntimeAttr? runtime_attr_merge_results
    }

    Array[File] vep_files = select_first([vep_vcf_files, vep_annotated_final_vcf])

    call step1.step1 as step1 {
        input:
            python_trio_sample_script=python_trio_sample_script,
            python_preprocess_script=python_preprocess_script,
            lcr_uri=lcr_uri,
            ped_sex_qc=ped_sex_qc,
            info_header=info_header,
            vep_files=vep_files,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix,
            hail_docker=hail_docker,
            hail_docker=hail_docker,
            qual_threshold=qual_threshold,
            sor_threshold_indel=sor_threshold_indel,
            sor_threshold_snv=sor_threshold_snv,
            readposranksum_threshold_indel=readposranksum_threshold_indel,
            readposranksum_threshold_snv=readposranksum_threshold_snv,
            qd_threshold_indel=qd_threshold_indel,
            qd_threshold_snv=qd_threshold_snv,
            mq_threshold=mq_threshold,
            exclude_gq_filters=exclude_gq_filters,
            merge_split_vcf=merge_split_vcf,
            shards_per_chunk=shards_per_chunk
    }

    call step2.step2 as step2 {
        input:
            merged_preprocessed_vcf_file=step1.merged_preprocessed_vcf_file,
            relatedness_qc=relatedness_qc,
            ped_sex_qc=ped_sex_qc,
            hail_docker=hail_docker
    }

    call step3.step3 as step3 {
        input:
            trio_uri=trio_uri,
            ped_sex_qc=ped_sex_qc,
            merged_preprocessed_vcf_file_filtered=step2.merged_preprocessed_vcf_file_filtered,
            hail_docker=hail_docker,
            cohort_prefix=cohort_prefix,
            trio_denovo_docker=trio_denovo_docker,
            uberSplit_v3_script=uberSplit_v3_script,
            batch_size=batch_size,
            subset_ped_script=subset_ped_script
    }
    call annotateHPandVAF.annotateHPandVAF as annotateHPandVAF {
        input:
            split_trio_vcfs=step3.split_trio_vcfs,
            vep_vcf_files=vep_files,
            hg38_reference=hg38_reference,
            hg38_reference_fai=hg38_reference_fai,
            hg38_reference_dict=hg38_reference_dict,
            jvarkit_docker=jvarkit_docker
    }

    call step4.step4 as step4 {
        input:
            ped_uri_trios=step3.ped_uri_trios,
            split_trio_vcfs=step3.split_trio_vcfs,
            get_sample_pedigree_script=get_sample_pedigree_script,
            trio_denovo_docker=trio_denovo_docker,
            minDQ=minDQ
    }

    call step5.step5 as step5 {
        input:
            ped_sex_qc=ped_sex_qc,
            split_trio_annot_vcfs=annotateHPandVAF.split_trio_annot_vcfs,
            trio_denovo_vcf=step4.trio_denovo_vcf,
            merge_vcf_to_tsv_fullQC_script=merge_vcf_to_tsv_fullQC_script,
            trio_denovo_docker=trio_denovo_docker,
            cohort_prefix=cohort_prefix
    }

    call annotateMPCandLOEUF.annotateMPCandLOEUF as annotateMPCandLOEUF {
        input:
            vcf_metrics_tsv=vcf_metrics_tsv,
            mpc_dir=mpc_dir,
            mpc_chr22_file=mpc_chr22_file,
            loeuf_file=loeuf_file,
            annotate_mpc_loeuf_script=annotate_mpc_loeuf_script,
            hail_docker=hail_docker
    }

    call step6.step6 as step6 {
        input:
            vcf_metrics_tsv_annot=annotateMPCandLOEUF.vcf_metrics_tsv_annot,
            vep_vcf_file=vep_files[0],
            AF_threshold=AF_threshold,
            AC_threshold=AC_threshold,
            csq_af_threshold=csq_af_threshold,
            filter_final_tsv_script=filter_final_tsv_script,
            hail_docker=hail_docker,
            sample_column=sample_column,
            genome_build=genome_build
    }

    call step7.step7 as step7 {
        input:
        vep_vcf_files=vep_files,
        lcr_uri=lcr_uri,
        ped_sex_qc=ped_sex_qc,
        meta_uri=meta_uri,
        trio_uri=trio_uri,
        vcf_metrics_tsv_final=step6.vcf_metrics_tsv_final,
        hg38_reference=hg38_reference,
        hg38_reference_dict=hg38_reference_dict,
        hg38_reference_fai=hg38_reference_fai,
        filter_rare_inherited_python_script=filter_rare_inherited_python_script,
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
        shards_per_chunk=shards_per_chunk,
        runtime_attr_merge_chunk=runtime_attr_merge_chunk,
        runtime_attr_filter_vcf=runtime_attr_filter_vcf,
        runtime_attr_merge_results=runtime_attr_merge_results,
        repetitive_regions_bed=repetitive_regions_bed,
        bagging_pu_source_script=bagging_pu_source_script,
        bagging_pu_rf_len_script=bagging_pu_rf_len_script,
        tsv_to_bed_script=tsv_to_bed_script,
        cohort_prefix=cohort_prefix,
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
        File meta_uri = step1.meta_uri
        File trio_uri = step1.trio_uri
        File merged_preprocessed_vcf_file = step1.merged_preprocessed_vcf_file
        File merged_preprocessed_vcf_idx = step1.merged_preprocessed_vcf_idx
        File merged_preprocessed_vcf_file_filtered = step2.merged_preprocessed_vcf_file_filtered
        File merged_preprocessed_sample_qc = step2.merged_preprocessed_sample_qc
        File ped_uri_trios = step3.ped_uri_trios
        Array[File] split_trio_vcfs = step3.split_trio_vcfs
        Array[File] split_trio_annot_vcfs = annotateHPandVAF.split_trio_annot_vcfs
        Array[File] trio_denovo_vcf = step4.trio_denovo_vcf
        File vcf_metrics_tsv = step5.vcf_metrics_tsv
        File vcf_metrics_tsv_annot = annotateMPCandLOEUF.vcf_metrics_tsv_annot
        File vcf_metrics_tsv_final = step6.vcf_metrics_tsv_final
        File vcf_metrics_tsv_final_pu = step7.vcf_metrics_tsv_final_pu
    }    
}