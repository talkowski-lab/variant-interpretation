
version 1.0

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

workflow step4 {
    input {
        Array[File] de_novo_results_sharded 
        Array[File] de_novo_vep_sharded 
        String sample_column
        String cohort_prefix
        String prioritize_csq_script
        String hail_docker
        RuntimeAttr? runtime_attr_merge_results
        RuntimeAttr? runtime_attr_prioritize
    }

    call helpers.mergeResultsPython as mergeDenovoResults {
        input:
            tsvs=de_novo_results_sharded,
            hail_docker=hail_docker,
            merged_filename=cohort_prefix + "_wes_final_denovo.tsv",
            input_size=size(de_novo_results_sharded, 'GB'),
            runtime_attr_override=runtime_attr_merge_results
    }

    call helpers.mergeResultsPython as mergeDenovoVEP {
        input:
            tsvs=de_novo_vep_sharded,
            hail_docker=hail_docker,
            merged_filename=cohort_prefix + "_wes_final_denovo_vep.tsv",
            input_size=size(de_novo_vep_sharded, 'GB'),
            runtime_attr_override=runtime_attr_merge_results
    }

    call prioritizeCSQ.mergeVEPIntoResults as mergeVEPIntoResults {
        input:
        de_novo_results=mergeDenovoResults.merged_tsv,
        de_novo_vep=mergeDenovoVEP.merged_tsv,
        cohort_prefix=cohort_prefix,
        hail_docker=hail_docker,
        runtime_attr_override=runtime_attr_merge_results
    }

    call prioritizeCSQ_og.annotateMostSevereCSQ as annotateMostSevereCSQ {
        input:
        vcf_metrics_tsv=mergeVEPIntoResults.de_novo_merged,
        prioritize_csq_script=prioritize_csq_script,
        hail_docker=hail_docker,
        sample_column=sample_column,
        runtime_attr_override=runtime_attr_prioritize
    }

    output {
        File de_novo_results = mergeDenovoResults.merged_tsv
        File de_novo_vep = mergeDenovoVEP.merged_tsv
        # prioritized CSQ
        File de_novo_merged = annotateMostSevereCSQ.vcf_metrics_tsv_prior_csq
    }
}

