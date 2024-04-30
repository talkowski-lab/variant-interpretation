version 1.0

import "wes-denovo-step-01-remote-sharded.wdl" as step1
import "wes-denovo-step-02-remote-sharded.wdl" as step2
import "wes-denovo-step-03-remote-sharded.wdl" as step3
import "wes-denovo-step-04-remote-sharded.wdl" as step4
import "wes-denovo-step-05-remote-sharded.wdl" as step5
import "wes-denovo-helpers.wdl" as helpers
import "wes-prioritize-csq.wdl" as prioritizeCSQ
import "prioritizeCSQ.wdl" as prioritizeCSQ_og

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
        Array[String] mt_uris
        File lcr_uri
        File ped_sex_qc
        File purcell5k
        File mpc_chr22_file
        File loeuf_file
        File eval_regions

        Boolean hail_autoscale
        String sample_column
        String bucket_id
        String mpc_dir
        String gnomad_ht_uri
        String cohort_prefix

        String hail_annotation_script
        String hail_basic_filtering_script
        String hail_denovo_filtering_script
        String prioritize_csq_script
        String final_filtering_script

        String hail_docker
        String sv_base_mini_docker
        String hail_docker

        Float min_child_ab=0.25
        Float min_dp_ratio=0.1
        Float min_gq=25
        Int vqslod_cutoff_snv=-20
        Int vqslod_cutoff_indel=-2
        Float af_threshold=0.005
        Float call_rate_threshold=0.8
        RuntimeAttr? runtime_attr_merge_results
        RuntimeAttr? runtime_attr_prioritize
    }

    scatter (mt_uri in mt_uris) {
        call helpers.getHailMTSize as getInputMTSize {
            input:
                mt_uri=mt_uri,
                hail_docker=hail_docker
        }
        call step1.hailAnnotateRemote as step1 {
            input:
                mt_uri=mt_uri,
                input_size=getInputMTSize.mt_size,
                ped_sex_qc=ped_sex_qc,
                purcell5k=purcell5k,
                mpc_chr22_file=mpc_chr22_file,
                mpc_dir=mpc_dir,
                gnomad_ht_uri=gnomad_ht_uri,
                bucket_id=bucket_id,
                cohort_prefix=cohort_prefix,
                hail_annotation_script=hail_annotation_script,
                hail_docker=hail_docker,
                hail_autoscale=hail_autoscale
        }

        call helpers.getHailMTSize as getStep1MTSize {
            input:
                mt_uri=step1.annot_mt,
                hail_docker=hail_docker
        }
        call step2.hailBasicFilteringRemote as step2 {
            input:
                lcr_uri=lcr_uri,
                annot_mt=step1.annot_mt,
                input_size=getStep1MTSize.mt_size,
                ped_sex_qc=ped_sex_qc,
                bucket_id=bucket_id,
                cohort_prefix=cohort_prefix,
                hail_basic_filtering_script=hail_basic_filtering_script,
                call_rate_threshold=call_rate_threshold,
                hail_docker=hail_docker
        }

        call helpers.getHailMTSize as getStep2MTSize {
            input:
                mt_uri=step2.filtered_mt,
                hail_docker=hail_docker
        }
    
        call step3.hailDenovoFilteringRemote as step3 {
            input:
                filtered_mt=step2.filtered_mt,
                input_size=getStep2MTSize.mt_size,
                ped_sex_qc=ped_sex_qc,
                bucket_id=bucket_id,
                cohort_prefix=cohort_prefix,
                loeuf_file=loeuf_file,
                hail_denovo_filtering_script=hail_denovo_filtering_script,
                hail_docker=hail_docker,
                min_child_ab=min_child_ab,
                min_dp_ratio=min_dp_ratio,
                min_gq=min_gq
        }
    }

    call step4.step4 as step4 {
        input:
        de_novo_results_sharded=step3.de_novo_results, 
        de_novo_vep_sharded=step3.de_novo_vep,
        sample_column=sample_column,
        cohort_prefix=cohort_prefix,
        prioritize_csq_script=prioritize_csq_script,
        hail_docker=hail_docker,
        runtime_attr_merge_results=runtime_attr_merge_results,
        runtime_attr_prioritize=runtime_attr_prioritize
    }

    call step5.step5 as step5 {
        input:
        de_novo_merged=step4.de_novo_merged,
        eval_regions=eval_regions,
        cohort_prefix=cohort_prefix,
        final_filtering_script=final_filtering_script,
        hail_docker=hail_docker,
        vqslod_cutoff_snv=vqslod_cutoff_snv,
        vqslod_cutoff_indel=vqslod_cutoff_indel,
        af_threshold=af_threshold
    }

    output {
        # step 1 output
        Array[String] annot_mt = step1.annot_mt
        Array[File] sample_qc_info = step1.sample_qc_info
        # step 2 output
        Array[String] filtered_mt = step2.filtered_mt
        Array[File] post_filter_sample_qc_info = step2.post_filter_sample_qc_info
        # step 3 output
        Array[String] de_novo_ht = step3.de_novo_ht
        Array[String] tdt_mt = step3.tdt_mt
        Array[String] tdt_parent_aware_mt = step3.tdt_parent_aware_mt
        # step4 output
        File de_novo_results = step4.de_novo_results
        File de_novo_vep = step4.de_novo_vep
        File de_novo_merged = step4.de_novo_merged
        # step5 output
        File de_novo_final = step5.de_novo_final
    }
}