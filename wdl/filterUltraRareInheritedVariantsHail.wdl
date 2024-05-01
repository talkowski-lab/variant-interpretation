version 1.0

import "mergeSplitVCF.wdl" as mergeSplitVCF
import "mergeVCFs.wdl" as mergeVCFs
import "wes-denovo-helpers.wdl" as helpers
import "downsampleVariantsfromTSV.wdl" as downsampleVariantsfromTSV

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow filterUltraRareInheritedVariantsHail {
    input {
        Array[File] vep_vcf_files
        File lcr_uri
        File ped_uri_trios
        File meta_uri
        File trio_uri
        File vcf_metrics_tsv_final
        File hg38_reference
        File hg38_reference_dict
        File hg38_reference_fai
        String python_trio_sample_script
        String filter_rare_inherited_python_script
        String jvarkit_docker
        String hail_docker
        String vep_hail_docker
        String sv_base_mini_docker
        String cohort_prefix
        Float AF_threshold=0.005
        Int AC_threshold=2
        Float csq_af_threshold=0.01
        Int gq_het_threshold=99
        Int gq_hom_ref_threshold=30
        Int qual_threshold=150
        Int shards_per_chunk=10
        RuntimeAttr? runtime_attr_merge_chunk
        RuntimeAttr? runtime_attr_filter_vcf
        RuntimeAttr? runtime_attr_merge_results
    }  

    # if (!defined(meta_uri)) {
    #     call makeTrioSampleFiles {
    #         input:
    #             python_trio_sample_script=python_trio_sample_script,
    #             ped_uri_trios=ped_uri_trios,
    #             cohort_prefix=cohort_prefix,
    #             hail_docker=hail_docker
    #     }        
    # }
    # File meta_uri_ = select_first([meta_uri, makeTrioSampleFiles.meta_uri])
    # File trio_uri_ = select_first([trio_uri, makeTrioSampleFiles.trio_uri])
    
    scatter (vcf_file in vep_vcf_files) {
        String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
        call filterUltraRareInheritedVariants as filterUltraRareInheritedVariants_sharded {
            input:
                vcf_file=vcf_file,
                lcr_uri=lcr_uri,
                ped_uri_trios=ped_uri_trios,
                meta_uri=meta_uri,
                trio_uri=trio_uri,
                filter_rare_inherited_python_script=filter_rare_inherited_python_script,
                vep_hail_docker=vep_hail_docker,
                cohort_prefix=basename(vcf_file, file_ext),
                AC_threshold=AC_threshold,
                AF_threshold=AF_threshold,
                csq_af_threshold=csq_af_threshold,
                gq_het_threshold=gq_het_threshold,
                gq_hom_ref_threshold=gq_hom_ref_threshold,
                qual_threshold=qual_threshold,
                runtime_attr_override=runtime_attr_filter_vcf
                }
    }
    call helpers.mergeResultsPython as mergeResults_sharded {
        input:
            tsvs=filterUltraRareInheritedVariants_sharded.ultra_rare_inherited_tsv,
            hail_docker=hail_docker,
            input_size=size(filterUltraRareInheritedVariants_sharded.ultra_rare_inherited_tsv, 'GB'),
            merged_filename=cohort_prefix+'_ultra_rare_variants.tsv.gz',
            runtime_attr_override=runtime_attr_merge_results
    }

    call downsampleVariantsfromTSV.downsampleVariantsfromTSV as downsampleVariantsfromTSV {
        input:
        reference_tsv=vcf_metrics_tsv_final,
        full_input_tsv=mergeResults_sharded.merged_tsv,
        hg38_reference=hg38_reference,
        hg38_reference_dict=hg38_reference_dict,
        hg38_reference_fai=hg38_reference_fai,
        jvarkit_docker=jvarkit_docker,
        vep_hail_docker=vep_hail_docker
    }

    output {
        File ultra_rare_inherited_tsv = mergeResults_sharded.merged_tsv
        File downsampled_ultra_rare_inherited_SNV = downsampleVariantsfromTSV.downsampled_tsv_SNV
        File downsampled_ultra_rare_inherited_Indel = downsampleVariantsfromTSV.downsampled_tsv_Indel
    }
}

task makeTrioSampleFiles {
    input {
        String python_trio_sample_script
        File ped_uri_trios
        String cohort_prefix
        String hail_docker
    }

    runtime {
        docker: hail_docker
    }

    command <<<
    curl ~{python_trio_sample_script} > python_trio_sample_script.py
    python3 python_trio_sample_script.py ~{ped_uri_trios} ~{cohort_prefix} 
    >>>
    
    output {
        File meta_uri = "~{cohort_prefix}_sample_list.txt"
        File trio_uri = "~{cohort_prefix}_trio_list.txt"
        File ped_uri_no_header = "~{cohort_prefix}_no_header.ped"
    }
}

task filterUltraRareInheritedVariants {
    input {
        File vcf_file
        File lcr_uri
        File ped_uri_trios
        File meta_uri
        File trio_uri
        String filter_rare_inherited_python_script
        String vep_hail_docker
        String cohort_prefix
        Int AC_threshold
        Float AF_threshold
        Float csq_af_threshold
        Int gq_het_threshold
        Int gq_hom_ref_threshold
        Int qual_threshold
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])

    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: vep_hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }


    command <<<
        curl ~{filter_rare_inherited_python_script} > filter_rare_variants.py
        python3.9 filter_rare_variants.py ~{lcr_uri} ~{ped_uri_trios} ~{meta_uri} ~{trio_uri} ~{vcf_file} \
        ~{cohort_prefix} ~{cpu_cores} ~{memory} ~{AC_threshold} ~{AF_threshold} ~{csq_af_threshold} \
        ~{gq_het_threshold} ~{gq_hom_ref_threshold} ~{qual_threshold} > stdout

        cp $(ls . | grep hail*.log) hail_log.txt
    >>>

    output {
        File hail_log = "hail_log.txt"
        File ultra_rare_inherited_tsv = cohort_prefix + '_ultra_rare_variants.tsv.gz'
    }
}
