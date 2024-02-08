version 1.0

import "compressHailMT.wdl" as compressHailMT
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
        String? mt_uri
        Float? mt_size
        File? annot_mt
        File ped_uri
        File purcell5k
        File mpc_chr22_file
        File loeuf_file
        Boolean sort_after_merge=false
        String mpc_dir
        String gnomad_ht_uri
        String cohort_prefix
        String hail_annotation_script
        String hail_basic_filtering_script
        String hail_denovo_filtering_script
        String hail_docker
        String sv_base_mini_docker
    }

    if (defined(mt_uri)) {
        call compressHailMT.compressMT as compressMT {
            input:
                mt_uri=select_first([mt_uri]),
                mt_size=select_first([mt_size]),
                hail_docker=hail_docker
        }
    }

    File annot_mt_ = select_first([compressMT.compressed_mt, annot_mt])
    
    call step2.hailBasicFiltering as step2 {
        input:
            annot_mt=annot_mt_,
            ped_uri=ped_uri,
            cohort_prefix=cohort_prefix,
            hail_basic_filtering_script=hail_basic_filtering_script,
            hail_docker=hail_docker
    }

    call step3.hailDenovoFiltering as step3 {
        input:
            filtered_mt=step2.filtered_mt,
            ped_uri=ped_uri,
            cohort_prefix=cohort_prefix,
            loeuf_file=loeuf_file,
            hail_denovo_filtering_script=hail_denovo_filtering_script,
            hail_docker=hail_docker
    }

    output {
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