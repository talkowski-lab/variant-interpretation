version 1.0

import "mergeSplitVCF.wdl" as mergeSplitVCF
import "mergeVCFs.wdl" as mergeVCFs

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow step2 {
    input {
        File merged_preprocessed_vcf_file
        File relatedness_qc
        File ped_sex_qc
        String hail_docker
        String remove_outliers_script
        RuntimeAttr? runtime_attr_remove_outliers
    }

    call removeOutliers {
        input:
        merged_preprocessed_vcf_file=merged_preprocessed_vcf_file,
        relatedness_qc=relatedness_qc,
        ped_sex_qc=ped_sex_qc,
        hail_docker=hail_docker,
        remove_outliers_script=remove_outliers_script,
        runtime_attr_override=runtime_attr_remove_outliers
    }

    output {
        File merged_preprocessed_vcf_file_filtered = removeOutliers.merged_preprocessed_vcf_file_filtered
    }
}

task removeOutliers {
    input {
        File merged_preprocessed_vcf_file
        File relatedness_qc
        File ped_sex_qc
        String hail_docker
        String remove_outliers_script
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(merged_preprocessed_vcf_file, "GB")
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

    command {
        curl ~{remove_outliers_script} > remove_outliers.py
        python3 remove_outliers.py ~{merged_preprocessed_vcf_file} ~{relatedness_qc} ~{ped_sex_qc} > stdout
    }

    output {
        File merged_preprocessed_vcf_file_filtered = basename(merged_preprocessed_vcf_file, '.vcf.gz') + '_removed_outliers.vcf.bgz' 
    }
}