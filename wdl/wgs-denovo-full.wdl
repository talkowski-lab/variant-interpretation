version 1.0

import "wgs-denovo-step-01.wdl" as step1
import "wgs-denovo-step-02.wdl" as step2
import "wgs-denovo-step-03.wdl" as step3
import "wgs-denovo-step-04.wdl" as step4
import "wgs-denovo-step-05.wdl" as step5
import "wgs-denovo-step-06.wdl" as step6
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
            Array[File]? vep_vcf_files
            Array[File]? vep_annotated_final_vcf
            String cohort_prefix
            String cohort_id
            String sv_base_mini_docker
            String trio_denovo_docker
            String hail_docker
            String vep_hail_docker
            String jvarkit_docker
            String sample_column
            Boolean exclude_gq_filters=false
            Boolean bad_header=false
            Boolean merge_split_vcf=false
            Int shards_per_chunk=10
            Int batch_size
            Float minDQ
            Float AF_threshold=0.005
            Int AC_threshold=2
            Float csq_af_threshold=0.01
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
            vep_hail_docker=vep_hail_docker,
            exclude_gq_filters=exclude_gq_filters,
            bad_header=bad_header,
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
            vep_annotated_final_vcf_single=vep_files[0],
            hg38_reference=hg38_reference,
            hg38_reference_fai=hg38_reference_fai,
            hg38_reference_dict=hg38_reference_dict,
            jvarkit_docker=jvarkit_docker
    }

    call step4.step4 as step4 {
        input:
            ped_sex_qc=step1.ped_sex_qc_no_header,
            split_trio_vcfs=step3.split_trio_vcfs,
            get_sample_pedigree_script=get_sample_pedigree_script,
            trio_denovo_docker=trio_denovo_docker,
            minDQ=minDQ
    }

    call step5.step5 as step5 {
        input:
            ped_sex_qc=ped_sex_qc,
            split_trio_vcfs=annotateHPandVAF.split_trio_annot_vcfs,
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
            AF_threshold=AF_threshold,
            AC_threshold=AC_threshold,
            csq_af_threshold=csq_af_threshold,
            filter_final_tsv_script=filter_final_tsv_script,
            hail_docker=hail_docker,
            sample_column=sample_column
    }

    output {
        File meta_uri = step1.meta_uri
        File trio_uri = step1.trio_uri
        File ped_sex_qc_no_header = step1.ped_sex_qc_no_header
        File merged_preprocessed_vcf_file = step1.merged_preprocessed_vcf_file
        File merged_preprocessed_vcf_idx = step1.merged_preprocessed_vcf_idx
        File merged_preprocessed_vcf_file_filtered = step2.merged_preprocessed_vcf_file_filtered
        File merged_preprocessed_sample_qc = step2.merged_preprocessed_sample_qc
        Array[File] split_trio_vcfs = step3.split_trio_vcfs
        Array[File] split_trio_annot_vcfs = annotateHPandVAF.split_trio_annot_vcfs
        Array[File] trio_denovo_vcf = step4.trio_denovo_vcf
        File vcf_metrics_tsv = step5.vcf_metrics_tsv
        File vcf_metrics_tsv_annot = annotateMPCandLOEUF.vcf_metrics_tsv_annot
        File vcf_metrics_tsv_final = step6.vcf_metrics_tsv_final
    }    
}