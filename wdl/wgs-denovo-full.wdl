version 1.0

import "wgs-denovo-step-01.wdl" as step1and2
import "wgs-denovo-step-03.wdl" as step3
import "wgs-denovo-step-04.wdl" as step4
import "wgs-denovo-step-05.wdl" as step5
import "annotateHPandVAF.wdl" as annotateHPandVAF

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
            File python_trio_sample_script
            File python_preprocess_script
            File subset_ped_python_script
            File uberSplit_v3_py
            File merge_vcf_to_tsv_fullQC_py
            File get_sample_pedigree_py
            File lcr_uri
            File ped_uri
            File hg38_reference
            File hg38_reference_fai
            File hg38_reference_dict
            File info_header
            Array[File]? vep_vcf_files
            Array[Array[File]]? vep_annotated_final_vcf
            String bucket_id
            String cohort_prefix
            String cohort_id
            String sv_base_mini_docker
            String trio_denovo_docker
            String hail_docker
            String vep_hail_docker
            String jvarkit_docker
            Boolean bad_header=false
            Int shards_per_chunk=0
            Int batch_size
            Float minDQ
    }

    if (defined(vep_annotated_final_vcf)) {
        Array[File] vep_annotated_final_vcf_arr = flatten(select_first([vep_annotated_final_vcf]))
    }
    Array[File] vep_files = select_first([vep_vcf_files, vep_annotated_final_vcf_arr])

    call step1and2.step1 as step1and2 {
        input:
            python_trio_sample_script=python_trio_sample_script,
            python_preprocess_script=python_preprocess_script,
            lcr_uri=lcr_uri,
            ped_uri=ped_uri,
            info_header=info_header,
            vep_files=vep_files,
            sv_base_mini_docker=sv_base_mini_docker,
            bucket_id=bucket_id,
            cohort_prefix=cohort_prefix,
            hail_docker=hail_docker,
            vep_hail_docker=vep_hail_docker,
            bad_header=bad_header,
            shards_per_chunk=shards_per_chunk
    }

    call step3.step3 as step3 {
        input:
            trio_uri=trio_uri,
            ped_uri=ped_uri,
            merged_preprocessed_vcf_file=step1and2.merged_preprocessed_vcf_file,
            hail_docker=hail_docker,
            cohort_prefix=cohort_prefix,
            trio_denovo_docker=trio_denovo_docker,
            uberSplit_v3_py=uberSplit_v3_py,
            batch_size=batch_size,
            subset_ped_python_script=subset_ped_python_script
    }
    call annotateHPandVAF.annotateHPandVAF as annotateHPandVAF {
        input:
            split_trio_vcfs=step3.split_trio_vcfs,
            vep_annotated_final_vcf_single=vep_files[0],
            hg38_reference=hg38_reference,
            hg38_reference_fai=hg38_reference_fai,
            hg38_reference_dict=hg38_reference_dict,
            jvarkit_docker=jvarkit_docker,
    }

    call step4.step4 as step4 {
        input:
            ped_uri=step1and2.ped_uri_no_header,
            split_trio_vcfs=step3.split_trio_vcfs,
            get_sample_pedigree_py=get_sample_pedigree_py,
            trio_denovo_docker=trio_denovo_docker,
            minDQ=minDQ
    }

    call step5.step5 as step5 {
        input:
            ped_uri=ped_uri,
            split_trio_vcfs=annotateHPandVAF.split_trio_annot_vcfs,
            trio_denovo_vcf=step4.trio_denovo_vcf,
            merge_vcf_to_tsv_fullQC_py=merge_vcf_to_tsv_fullQC_py,
            trio_denovo_docker=trio_denovo_docker,
            cohort_prefix=cohort_prefix
    }

    output {
        File meta_uri = step1and2.meta_uri
        File trio_uri = step1and2.trio_uri
        File ped_uri_no_header = step1and2.ped_uri_no_header
        File merged_preprocessed_vcf_file = step1and2.merged_preprocessed_vcf_file
        File merged_preprocessed_vcf_idx = step1and2.merged_preprocessed_vcf_idx
        Array[File] split_trio_vcfs = step3.split_trio_vcfs
        Array[File] split_trio_annot_vcfs = annotateHPandVAF.split_trio_annot_vcfs
        Array[File] trio_denovo_vcf = step4.trio_denovo_vcf
        File vcf_metrics_tsv = step5.vcf_metrics_tsv
    }    
}