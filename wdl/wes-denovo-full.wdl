version 1.0

import "wes-denovo-step-01.wdl" as step1
import "wes-denovo-step-02.wdl" as step2
import "wes-denovo-step-03.wdl" as step3
import "mergeVCFs.wdl" as mergeVCFs

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow hailDenovoWES {
    input {
        File? vcf_file
        Array[File]? vcf_files
        File ped_uri
        File purcell5k
        File mpc_chr22_file
        File loeuf_file
        Boolean sort_after_merge=false
        String bucket_id=false
        String mpc_dir
        String gnomad_ht_uri
        String cohort_prefix
        String hail_annotation_script
        String hail_basic_filtering_script
        String hail_denovo_filtering_script
        String hail_docker
        String sv_base_mini_docker
    }

    if (defined(vcf_files)) {
        call mergeVCFs.mergeVCFs as mergeVCFs {
            input:
                vcf_files=select_first([vcf_files]),
                sv_base_mini_docker=sv_base_mini_docker,
                cohort_prefix=cohort_prefix,
                sort_after_merge=sort_after_merge
        }
    }

    File vcf_file_ = select_first([vcf_file, mergeVCFs.merged_vcf_file])

    call step1.hailAnnotate as step1 {
        input:
            vcf_file=vcf_file_,
            ped_uri=ped_uri,
            purcell5k=purcell5k,
            mpc_chr22_file=mpc_chr22_file,
            mpc_dir=mpc_dir,
            gnomad_ht_uri=gnomad_ht_uri,
            bucket_id=bucket_id,
            cohort_prefix=cohort_prefix,
            hail_annotation_script=hail_annotation_script,
            hail_docker=hail_docker
    }

    call step2.hailBasicFiltering as step2 {
        input:
            annot_mt=step1.annot_mt,
            ped_uri=ped_uri,
            cohort_prefix=cohort_prefix,
            bucket_id=bucket_id,
            hail_basic_filtering_script=hail_basic_filtering_script,
            hail_docker=hail_docker
    }

    call step3.hailDenovoFiltering as step3 {
        input:
            filtered_mt=step2.filtered_mt,
            ped_uri=ped_uri,
            cohort_prefix=cohort_prefix,
            bucket_id=bucket_id,
            loeuf_file=loeuf_file,
            hail_denovo_filtering_script=hail_denovo_filtering_script,
            hail_docker=hail_docker
    }

    output {
        # step 1 output
        File annot_mt = step1.annot_mt
        File sample_qc_info = step1.sample_qc_info
        File pca_score_table_5k = step1.pca_score_table_5k
        File pca_loading_table_5k = step1.pca_loading_table_5k
        # step 2 output
        File filtered_mt = step2.filtered_mt
        File post_filter_sample_qc_info = step2.post_filter_sample_qc_info
        # step 3 output
        File de_novo_results = step3.de_novo_results
        File de_novo_vep = step3.de_novo_vep
        File de_novo_ht = step3.de_novo_ht
        File tdt_mt = step3.tdt_mt
        File tdt_parent_aware_mt = step3.tdt_parent_aware_mt
    }
}
