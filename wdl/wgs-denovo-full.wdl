version 1.0

import "wgs-denovo-step-01.wdl" as step1and2
import "wgs-denovo-step-03.wdl" as step3
import "wgs-denovo-step-04.wdl" as step4
import "wgs-denovo-step-05.wdl" as step5

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
            # file can be a list of vcf files or just one vcf file
            File file
            File python_trio_sample_script
            File python_preprocess_script
            File uberSplit_v3_py
            File merge_vcf_to_tsv_fullQC_py
            File lcr_uri
            File ped_uri
            File ped_uri_no_header
            Array[Array[File]] vcf_uri_list
            String bucket_id
            String cohort_prefix
            String sv_base_mini_docker
            String trio_denovo_docker
            Int batch_size
    }

    call step1and2.step1 as step1and2 {
        input:
            file=file,
            python_trio_sample_script=python_trio_sample_script,
            python_preprocess_script=python_preprocess_script,
            lcr_uri=lcr_uri,
            ped_uri=ped_uri,
            vcf_uri_list=vcf_uri_list,
            sv_base_mini_docker=sv_base_mini_docker,
            bucket_id=bucket_id,
            cohort_prefix=cohort_prefix,
            trio_denovo_docker=trio_denovo_docker
    }

    call step3.step3 as step3 {
        input:
            ped_uri=ped_uri,
            merged_preprocessed_vcf_files=step1and2.merged_preprocessed_vcf_files,
            trio_denovo_docker=trio_denovo_docker,
            uberSplit_v3_py=uberSplit_v3_py,
            batch_size=batch_size
    }

    call step4.step4 as step4 {
        input:
            ped_uri=ped_uri_no_header,
            split_trio_vcfs=step3.split_trio_vcfs,
            trio_denovo_docker=trio_denovo_docker
    }

    call step5.step5 as step5 {
        input:
            ped_uri=ped_uri,
            split_trio_vcfs=step3.split_trio_vcfs,
            trio_denovo_vcf=step4.trio_denovo_vcf,
            merge_vcf_to_tsv_fullQC_py=merge_vcf_to_tsv_fullQC_py,
            trio_denovo_docker=trio_denovo_docker,
            cohort_prefix=cohort_prefix
    }

    output {
        Array[File] merged_preprocessed_vcf_files = step1and2.merged_preprocessed_vcf_files
        Array[File] merged_preprocessed_vcf_idx = step1and2.merged_preprocessed_vcf_idx
        Array[Array[File]] split_trio_vcfs = step3.split_trio_vcfs
        Array[File] stats_files = step3.stats_files
        Array[Array[File]] trio_denovo_vcf = step4.trio_denovo_vcf
        File vcf_metrics_tsv = step5.vcf_metrics_tsv
    }    
}