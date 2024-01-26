version 1.0

import "wes-denovo-step-01.wdl" as step1
import "wes-denovo-step-02.wdl" as step2
import "wes-denovo-step-03.wdl" as step3

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
        File vcf_file
        File corrected_ped
        File purcell5k
        File mpc_chr22_file
        File loeuf_file
        String mpc_dir
        String gnomad_ht_uri
        String cohort_prefix
        String hail_annotation_script
        String hail_basic_filtering_script
        String hail_denovo_filtering_script
        String hail_docker
    }

    call step1.hailAnnotate as step1 {
        input:
            vcf_file=vcf_file,
            corrected_ped=corrected_ped,
            purcell5k=purcell5k,
            mpc_chr22_file=mpc_chr22_file,
            mpc_dir=mpc_dir,
            gnomad_ht_uri=gnomad_ht_uri,
            cohort_prefix=cohort_prefix,
            hail_annotation_script=hail_annotation_script,
            hail_docker=hail_docker
    }

    call step2.hailBasicFiltering as step2 {
        input:
            annot_mt=step1.annot_mt,
            corrected_ped=corrected_ped,
            cohort_prefix=cohort_prefix,
            hail_basic_filtering_script=hail_basic_filtering_script,
            hail_docker=hail_docker
    }

    call step3.hailDenovoFiltering as step3 {
        input:
            filtered_mt=step2.filtered_mt,
            corrected_ped=corrected_ped,
            cohort_prefix=cohort_prefix,
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