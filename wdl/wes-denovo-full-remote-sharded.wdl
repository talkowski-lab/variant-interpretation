version 1.0

import "wes-denovo-step-01.wdl" as step1
import "wes-denovo-step-02.wdl" as step2
import "wes-denovo-step-03.wdl" as step3
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
        File ped_uri
        File purcell5k
        File mpc_chr22_file
        File loeuf_file
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
        String hail_docker
        String sv_base_mini_docker
        Float call_rate_threshold=0.8
    }

    scatter (mt_uri in mt_uris) {
        call helpers.getHailMTSize as getInputMTSize {
            input:
                mt_uri=mt_uri,
                hail_docker=hail_docker
        }
        call step1.hailAnnotateRemote as step1 {
            input:
                vcf_file=mt_uri,
                input_size=getInputMTSize.mt_size,
                ped_uri=ped_uri,
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
                ped_uri=ped_uri,
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
                ped_uri=ped_uri,
                bucket_id=bucket_id,
                cohort_prefix=cohort_prefix,
                loeuf_file=loeuf_file,
                hail_denovo_filtering_script=hail_denovo_filtering_script,
                hail_docker=hail_docker
        }
    }

    call helpers.mergeResultsPython as mergeDenovoResults {
        input:
            tsvs=step3.de_novo_results,
            hail_docker=hail_docker,
            merged_filename=cohort_prefix + "_wes_final_denovo.tsv",
            input_size=size(step3.de_novo_results, 'GB')
    }

    call helpers.mergeResultsPython as mergeDenovoVEP {
        input:
            tsvs=step3.de_novo_vep,
            hail_docker=hail_docker,
            merged_filename=cohort_prefix + "_wes_final_denovo_vep.tsv",
            input_size=size(step3.de_novo_vep, 'GB')
    }

    call prioritizeCSQ.mergeVEPIntoResults as mergeVEPIntoResults {
        input:
        de_novo_results=mergeDenovoResults.merged_tsv,
        de_novo_vep=mergeDenovoVEP.merged_tsv,
        cohort_prefix=cohort_prefix,
        hail_docker=hail_docker
    }

    call prioritizeCSQ_og.annotateMostSevereCSQ as annotateMostSevereCSQ {
        input:
        vcf_metrics_tsv=mergeVEPIntoResults.de_novo_merged,
        prioritize_csq_script=prioritize_csq_script,
        hail_docker=hail_docker,
        sample_column=sample_column
    }

    output {
        # step 1 output
        Array[String] annot_mt = step1.annot_mt
        Array[File] sample_qc_info = step1.sample_qc_info
        # step 2 output
        Array[String] filtered_mt = step2.filtered_mt
        Array[File] post_filter_sample_qc_info = step2.post_filter_sample_qc_info
        # step 3 output
        File de_novo_results = mergeDenovoResults.merged_tsv
        File de_novo_vep = mergeDenovoVEP.merged_tsv
        Array[String] de_novo_ht = step3.de_novo_ht
        Array[String] tdt_mt = step3.tdt_mt
        Array[String] tdt_parent_aware_mt = step3.tdt_parent_aware_mt
        # prioritized CSQ
        File de_novo_merged = annotateMostSevereCSQ.vcf_metrics_tsv_prior_csq
    }
}