version 1.0 

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow step5 {
    input {
        File ped_sex_qc
        Array[File] split_trio_annot_vcfs  # for input directory (from step 4)
        Array[File] trio_denovo_vcf  # for output directory
        String merge_vcf_to_tsv_fullQC_script
        String trio_denovo_docker
        String cohort_prefix
        RuntimeAttr? runtime_attr_step5
    }

    call merge_vcf_to_tsv_fullQC {
        input:
            ped_sex_qc=ped_sex_qc,
            split_trio_annot_vcfs=split_trio_annot_vcfs,
            trio_denovo_vcf=trio_denovo_vcf,
            merge_vcf_to_tsv_fullQC_script=merge_vcf_to_tsv_fullQC_script,
            trio_denovo_docker=trio_denovo_docker,
            cohort_prefix=cohort_prefix,
            runtime_attr_override=runtime_attr_step5
    }

    output {
        File vcf_metrics_tsv = merge_vcf_to_tsv_fullQC.output_tsv
    }
}

task merge_vcf_to_tsv_fullQC {
    input {
        File ped_sex_qc
        String merge_vcf_to_tsv_fullQC_script
        Array[File] split_trio_annot_vcfs 
        Array[File] trio_denovo_vcf
        String trio_denovo_docker
        String cohort_prefix
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(split_trio_annot_vcfs, "GB") + size(trio_denovo_vcf, "GB")
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
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: trio_denovo_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command {
        input_dir=$(dirname ~{split_trio_annot_vcfs[0]})
        output_dir=$(dirname ~{trio_denovo_vcf[0]})
        curl ~{merge_vcf_to_tsv_fullQC_script} > merge_vcf_to_tsv_fullQC.py
        python3 merge_vcf_to_tsv_fullQC.py -d $output_dir -i $input_dir -p ~{ped_sex_qc} -o ~{cohort_prefix}_dnm.tsv > stdout
    }

    output {
        File output_tsv = cohort_prefix + '_dnm.tsv'
    }
}