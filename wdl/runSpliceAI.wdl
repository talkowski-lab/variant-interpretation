version 1.0

import "mergeVCFs.wdl" as mergeVCFs

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? gpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow filterRunSpliceAI {
    input {
        Array[File] vep_vcf_files
        File gene_annotation_file
        File ref_fasta

        String cohort_prefix
        String filter_rare_splice_impact_script
        String hail_docker
        String sv_base_mini_docker
        String spliceAI_docker

        Int ac_threshold=20
        Float gnomad_af_threshold=0.05
        Int max_distance=50
        Boolean mask=false
        Boolean sort_after_merge=false
    }
    
    scatter (vcf_file in vep_vcf_files) {
        call filterRareSpliceImpactVariants {
            input:
            vcf_file=vcf_file,
            filter_rare_splice_impact_script=filter_rare_splice_impact_script,
            hail_docker=hail_docker,
            ac_threshold=ac_threshold,
            gnomad_af_threshold=gnomad_af_threshold
        }

        call runSpliceAI {
            input:
            vcf_file=filterRareSpliceImpactVariants.filtered_vcf,
            vcf_idx=filterRareSpliceImpactVariants.filtered_vcf_idx,
            ref_fasta=ref_fasta,
            gene_annotation_file=gene_annotation_file,
            spliceAI_docker=spliceAI_docker,
            max_distance=max_distance,
            mask=mask
        }
    }

    call mergeVCFs.mergeVCFs as mergeVCFs {
        input:  
            vcf_files=runSpliceAI.spliceAI_vcf,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix + '_filtered_spliceAI',
            sort_after_merge=sort_after_merge        
    }

    output {
        File filtered_spliceAI_vcf = mergeVCFs.merged_vcf_file
        File filtered_spliceAI_vcf_idx = mergeVCFs.merged_vcf_idx
    }
}
