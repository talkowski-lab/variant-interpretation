version 1.0

import "mergeSplitVCF.wdl" as mergeSplitVCF
import "mergeVCFs.wdl" as mergeVCFs
import "wes-denovo-helpers.wdl" as helpers
import "downsampleVariantsfromTSV.wdl" as downsampleVariantsfromTSV
import "prioritizeCSQ.wdl" as prioritizeCSQ

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow filterUltraRareParentsVariantsHail {
    input {
        Array[File] annot_vcf_files
        File lcr_uri
        File ped_sex_qc
        File meta_uri
        File trio_uri
        File vcf_metrics_tsv_final
        File remove_regions_bed
        File hg38_reference
        File hg38_reference_dict
        File hg38_reference_fai
        # String python_trio_sample_script
        String filter_rare_parents_python_script
        String jvarkit_docker
        String hail_docker
        String sv_base_mini_docker
        String cohort_prefix
        Boolean sort_after_merge=false
        Float AF_threshold=0.005
        Int AC_threshold=1
        Float csq_af_threshold=0.00001
        Int gq_het_threshold=99
        Int gq_hom_ref_threshold=30
        Int qual_threshold=150  # ~30 for DRAGEN
        Float sor_threshold_indel=3.0
        Float sor_threshold_snv=2.5
        Float readposranksum_threshold_indel=-1.7
        Float readposranksum_threshold_snv=-1.4
        Float qd_threshold_indel=4.0
        Float qd_threshold_snv=3.0
        Float mq_threshold=50
        Int shards_per_chunk=10
        String genome_build='GRCh38'

        # for prioritizeCSQ
        String prioritize_csq_script
        String sample_column='SAMPLE'

        # for downsampling
        Int chunk_size=100000
        Float snv_scale=1
        Float indel_scale=1
        Boolean prioritize_coding=true
        Boolean prioritize_gnomad=true

        RuntimeAttr? runtime_attr_filter_vcf
        RuntimeAttr? runtime_attr_merge_results
        RuntimeAttr? runtime_attr_prioritize
        RuntimeAttr? runtime_attr_downsample
    }  

    # if (!defined(meta_uri)) {
    #     call makeTrioSampleFiles {
    #         input:
    #             python_trio_sample_script=python_trio_sample_script,
    #             ped_sex_qc=ped_sex_qc,
    #             cohort_prefix=cohort_prefix,
    #             hail_docker=hail_docker
    #     }        
    # }
    # File meta_uri_ = select_first([meta_uri, makeTrioSampleFiles.meta_uri])
    # File trio_uri_ = select_first([trio_uri, makeTrioSampleFiles.trio_uri])

    scatter (vcf_file in annot_vcf_files) {
        String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
        call filterRareParentsVariants as filterRareParentsVariants_sharded {
            input:
                vcf_file=vcf_file,
                lcr_uri=lcr_uri,
                ped_sex_qc=ped_sex_qc,
                meta_uri=meta_uri,
                trio_uri=trio_uri,
                filter_rare_parents_python_script=filter_rare_parents_python_script,
                hail_docker=hail_docker,
                cohort_prefix=basename(vcf_file, file_ext),
                AC_threshold=AC_threshold,
                AF_threshold=AF_threshold,
                csq_af_threshold=csq_af_threshold,
                gq_het_threshold=gq_het_threshold,
                gq_hom_ref_threshold=gq_hom_ref_threshold,
                qual_threshold=qual_threshold,
                sor_threshold_indel=sor_threshold_indel,
                sor_threshold_snv=sor_threshold_snv,
                readposranksum_threshold_indel=readposranksum_threshold_indel,
                readposranksum_threshold_snv=readposranksum_threshold_snv,
                qd_threshold_indel=qd_threshold_indel,
                qd_threshold_snv=qd_threshold_snv,
                mq_threshold=mq_threshold,
                genome_build=genome_build,
                runtime_attr_override=runtime_attr_filter_vcf
                }
    }
    call helpers.mergeResultsPython as mergeResults_sharded {
        input:
            tsvs=filterRareParentsVariants_sharded.ultra_rare_parents_tsv,
            hail_docker=hail_docker,
            input_size=size(filterRareParentsVariants_sharded.ultra_rare_parents_tsv, 'GB'),
            merged_filename=cohort_prefix+'_ultra_rare_parents_variants.tsv.gz',
            runtime_attr_override=runtime_attr_merge_results
    }

    call prioritizeCSQ.annotateMostSevereCSQ as prioritizeCSQ {
        input:
        vcf_metrics_tsv=mergeResults_sharded.merged_tsv,
        vep_vcf_file=annot_vcf_files[0],
        hail_docker=hail_docker,
        prioritize_csq_script=prioritize_csq_script,
        sample_column=sample_column,
        genome_build=genome_build,
        runtime_attr_override=runtime_attr_prioritize
    }

    call downsampleVariantsfromTSV.downsampleVariantsfromTSV as downsampleVariantsfromTSV {
        input:
        reference_tsv=vcf_metrics_tsv_final,
        full_input_tsv=prioritizeCSQ.vcf_metrics_tsv_prior_csq,
        remove_regions_bed=remove_regions_bed,
        hg38_reference=hg38_reference,
        hg38_reference_dict=hg38_reference_dict,
        hg38_reference_fai=hg38_reference_fai,
        jvarkit_docker=jvarkit_docker,
        hail_docker=hail_docker,
        chunk_size=chunk_size,
        snv_scale=snv_scale,
        indel_scale=indel_scale,
        prioritize_gnomad=prioritize_gnomad,
        runtime_attr_downsample=runtime_attr_downsample
    }

    output {
        File ultra_rare_parents_tsv = prioritizeCSQ.vcf_metrics_tsv_prior_csq
        File downsampled_ultra_rare_parents_SNV = downsampleVariantsfromTSV.downsampled_tsv_SNV
        File downsampled_ultra_rare_parents_Indel = downsampleVariantsfromTSV.downsampled_tsv_Indel
    }
}

task makeTrioSampleFiles {
    input {
        String python_trio_sample_script
        File ped_sex_qc
        String cohort_prefix
        String hail_docker
    }

    runtime {
        docker: hail_docker
    }

    command <<<
    curl ~{python_trio_sample_script} > python_trio_sample_script.py
    python3 python_trio_sample_script.py ~{ped_sex_qc} ~{cohort_prefix} 
    >>>
    
    output {
        File meta_uri = "~{cohort_prefix}_sample_list.txt"
        File trio_uri = "~{cohort_prefix}_trio_list.txt"
        File ped_uri_no_header = "~{cohort_prefix}_no_header.ped"
    }
}

task filterRareParentsVariants {
    input {
        File vcf_file
        File lcr_uri
        File ped_sex_qc
        File meta_uri
        File trio_uri
        String filter_rare_parents_python_script
        String hail_docker
        String cohort_prefix
        Int AC_threshold
        Float AF_threshold
        Float csq_af_threshold
        Int gq_het_threshold
        Int gq_hom_ref_threshold
        Int qual_threshold
        Float sor_threshold_indel
        Float sor_threshold_snv
        Float readposranksum_threshold_indel
        Float readposranksum_threshold_snv
        Float qd_threshold_indel
        Float qd_threshold_snv
        Float mq_threshold
        String genome_build
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
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }


    command <<<
        set -eou pipefail
        curl ~{filter_rare_parents_python_script} > filter_rare_variants.py
        python3 filter_rare_variants.py ~{lcr_uri} ~{ped_sex_qc} ~{meta_uri} ~{trio_uri} ~{vcf_file} \
        ~{cohort_prefix} ~{cpu_cores} ~{memory} ~{AC_threshold} ~{AF_threshold} ~{csq_af_threshold} \
        ~{gq_het_threshold} ~{gq_hom_ref_threshold} ~{qual_threshold} ~{sor_threshold_indel} ~{sor_threshold_snv} \
        ~{readposranksum_threshold_indel} ~{readposranksum_threshold_snv} ~{qd_threshold_indel} ~{qd_threshold_snv} \
        ~{mq_threshold} ~{genome_build} > stdout

        cp $(ls . | grep hail*.log) hail_log.txt
    >>>

    output {
        File hail_log = "hail_log.txt"
        File ultra_rare_parents_tsv = cohort_prefix + '_ultra_rare_parents_variants.tsv.gz'
    }
}
