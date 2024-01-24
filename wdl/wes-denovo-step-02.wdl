version 1.0

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
        File corrected_ped
        String annot_mt
        String cohort_prefix
        String bucket_id
        String hail_basic_filtering_script
        String hail_docker
    }

    call hailBasicFiltering {
        input:
            annot_mt=annot_mt,
            corrected_ped=corrected_ped,
            cohort_prefix=cohort_prefix,
            bucket_id=bucket_id,
            hail_basic_filtering_script=hail_basic_filtering_script,
            hail_docker=hail_docker
    }

    output {
        File filtered_mt = hailBasicFiltering.filtered_mt
        File post_filter_sample_qc_info = hailBasicFiltering.post_filter_sample_qc_info
    }
}

task hailBasicFiltering {
    input {
        File corrected_ped
        String annot_mt
        String cohort_prefix
        String bucket_id
        String hail_basic_filtering_script
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(corrected_ped, "GB")
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
        curl ~{hail_basic_filtering_script} > hail_basic_filtering_script.py
        python3 hail_basic_filtering_script.py ~{annot_mt} ~{cohort_prefix} ~{corrected_ped} ~{bucket_id} \
        ~{cpu_cores} ~{memory}
    }

    output {
        String filtered_mt = "~{bucket_id}/hail/~{cohort_prefix}_wes_denovo_basic_filtering.mt"
        File post_filter_sample_qc_info = "~{cohort_prefix}_wes_final_annot_post_filter_qc_info.txt"
    }
}