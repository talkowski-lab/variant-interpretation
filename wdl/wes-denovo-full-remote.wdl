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
        String mt_uri
        Float mt_size
        File ped_uri
        File purcell5k
        File mpc_chr22_file
        File loeuf_file
        String bucket_id
        String mpc_dir
        String gnomad_ht_uri
        String cohort_prefix
        String hail_annotation_script
        String hail_basic_filtering_script
        String hail_denovo_filtering_script
        String hail_docker
        String sv_base_mini_docker
    }

    call step1.hailAnnotateRemote as step1 {
        input:
            mt_uri=mt_uri,
            input_size=mt_size,
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

    call step2.hailBasicFilteringRemote as step2 {
        input:
            annot_mt=mt_uri,
            input_size=mt_size,
            ped_uri=ped_uri,
            bucket_id=bucket_id,
            cohort_prefix=cohort_prefix,
            hail_basic_filtering_script=hail_basic_filtering_script,
            hail_docker=hail_docker
    }

    call step3.hailDenovoFilteringRemote as step3 {
        input:
            filtered_mt=step2.filtered_mt,
            input_size=mt_size,
            ped_uri=ped_uri,
            bucket_id=bucket_id,
            cohort_prefix=cohort_prefix,
            loeuf_file=loeuf_file,
            hail_denovo_filtering_script=hail_denovo_filtering_script,
            hail_docker=hail_docker
    }

    output {
        # step 2 output
        String filtered_mt = step2.filtered_mt
        File post_filter_sample_qc_info = step2.post_filter_sample_qc_info
        # step 3 output
        File de_novo_results = step3.de_novo_results
        File de_novo_vep = step3.de_novo_vep
        String de_novo_ht = step3.de_novo_ht
        String tdt_mt = step3.tdt_mt
        String tdt_parent_aware_mt = step3.tdt_parent_aware_mt
    }
}