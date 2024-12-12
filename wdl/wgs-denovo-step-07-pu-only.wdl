version 1.0

# assuming filterUltraRareInheritedVariantsHail and filterUltraRareParentsVariantsHail already run
import  "flagRepetitiveRegions.wdl" as flagRepetitiveRegions
import "wgs-denovo-bagging-pu-rf-len.wdl" as BaggingPU_RF
import "filterUltraRareInheritedVariantsHail.wdl" as filterUltraRareInheritedVariantsHail
import "filterUltraRareParentsVariantsHail.wdl" as filterUltraRareParentsVariantsHail 
import "wes-denovo-helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow step7 {
    input {
        File downsampled_ultra_rare_inherited
        File downsampled_ultra_rare_parents
        File lcr_uri
        File ped_sex_qc
        File meta_uri
        File trio_uri
        File vcf_metrics_tsv_final
        File hg38_reference
        File hg38_reference_dict
        File hg38_reference_fai
        String jvarkit_docker
        String hail_docker
        String sv_base_mini_docker
        String cohort_prefix
    
        # 12/12/2024 NEW
        Boolean batch_coding_only=true
        Int batch_size=10000
        RuntimeAttr? runtime_attr_batch
        RuntimeAttr? runtime_attr_merge_results

        File repetitive_regions_bed
        String var_type  # Indel or SNV
        String bagging_pu_source_script
        String bagging_pu_rf_len_script
        String tsv_to_bed_script
        String cohort_prefix
        String metric='fp_fn_ratio'
        Array[String] sample_features=["GQ_parent", "AB_sample", "DPC_sample", "DPC_parent", "PL_sample_0.0", "PL_sample_1.1"]
        Array[String] variant_features=["MQ", "FS", "BaseQRankSum", "SOR", "LEN", "ReadPosRankSum", "DP", "QD", "VQSLOD"]
        Float vqslod_cutoff=-10
        Int n_estimators_rf=100
        Int n_bag=10
        Int n_jobs=-1
        Boolean filter_pass_before=false
        RuntimeAttr? runtime_attr_bagging_pu
    }

    call splitIntoBatches {
        input:
            vcf_metrics_tsv_final=vcf_metrics_tsv_final,
            var_type=var_type,
            hail_docker=hail_docker,
            batch_coding_only=batch_coding_only,
            batch_size=batch_size,
            runtime_attr_override=runtime_attr_batch
    }

    scatter (batch_tsv in splitIntoBatches.batches) {
        call BaggingPU_RF.BaggingPU_RF as BaggingPU_RF {
            input:
            vcf_metrics_tsv_final=batch_tsv,
            ultra_rare_inherited_tsv=downsampled_ultra_rare_inherited,
            ultra_rare_parents_tsv=downsampled_ultra_rare_parents,
            repetitive_regions_bed=repetitive_regions_bed,
            var_type=var_type,
            bagging_pu_source_script=bagging_pu_source_script,
            bagging_pu_rf_len_script=bagging_pu_rf_len_script,
            tsv_to_bed_script=tsv_to_bed_script,
            cohort_prefix=cohort_prefix,
            sv_base_mini_docker=sv_base_mini_docker,
            hail_docker=hail_docker,
            metric=metric,
            sample_features=sample_features,
            variant_features=variant_features,
            vqslod_cutoff=vqslod_cutoff,
            n_estimators_rf=n_estimators_rf,
            n_bag=n_bag,
            n_jobs=n_jobs,
            filter_pass_before=filter_pass_before,
            runtime_attr_bagging_pu=runtime_attr_bagging_pu
        }
    }

    call helpers.mergeResultsPython as mergeBaggingPU {
        input:
            tsvs=BaggingPU_RF.vcf_metrics_tsv_final_pu,
            hail_docker=hail_docker,
            input_size=size(BaggingPU_RF.vcf_metrics_tsv_final_pu, 'GB'),
            merged_filename="~{basename(vcf_metrics_tsv_final, '.tsv.gz')}_pu_~{var_type}.tsv.gz",
            runtime_attr_override=runtime_attr_merge_results
    }

    output {
        File vcf_metrics_tsv_final_pu = mergeBaggingPU.merged_tsv
        Array[File] pu_feature_importances_plot = BaggingPU_RF.pu_feature_importances_plot
    }
}

task splitIntoBatches {
    input {
        File vcf_metrics_tsv_final
        String var_type
        String hail_docker

        Boolean batch_coding_only
        Int batch_size

        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_metrics_tsv_final, "GB")
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

    String file_ext = if sub(basename(vcf_metrics_tsv_final), '.tsv.gz', '')!=basename(vcf_metrics_tsv_final) then '.tsv.gz' else '.tsv'

    command <<<
    cat <<EOF > split_into_batches.py
    import pandas as pd
    import numpy as np
    import os
    import sys
    import ast

    file = sys.argv[1]
    batch_coding_only = ast.literal_eval(sys.argv[2].capitalize())
    batch_size = int(sys.argv[3])
    var_type = sys.argv[4]
    file_ext = sys.argv[5]

    df = pd.concat(pd.read_csv(file, sep='\t', chunksize=100_000))
    df = df[df.TYPE==var_type].copy()

    base_filename = os.path.basename(file).split(file_ext)[0]
    i=0
    if batch_coding_only:  # save all coding to one batch
        df[df.isCoding].to_csv(f"{base_filename}.shard_{i}.all_coding_only.{file_ext}", sep='\t', index=False)
        df = df[~df.isCoding].copy()
        i+=1

    while df.shape[0]>0:
        batch = df.sample(min(batch_size, df.shape[0]))
        batch.to_csv(f"{base_filename}.shard_{i}.{file_ext}", sep='\t', index=False)
        df = df.drop(batch.index)    
        i+=1    
    EOF
    python3 split_into_batches.py ~{vcf_metrics_tsv_final} ~{batch_coding_only} ~{batch_size} ~{var_type} ~{file_ext}
    >>>

    String base_filename = basename(vcf_metrics_tsv_final, file_ext)
    output {
        Array[File] batches = glob("~{base_filename}.shard_*")
    }
}