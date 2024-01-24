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
        String filtered_mt
        String cohort_prefix
        String bucket_id
        String hail_denovo_filtering_script
        String hail_docker
    }

    call hailDenovoFiltering {
        input:
            filtered_mt=filtered_mt,
            corrected_ped=corrected_ped,
            cohort_prefix=cohort_prefix,
            bucket_id=bucket_id,
            hail_denovo_filtering_script=hail_denovo_filtering_script,
            hail_docker=hail_docker
    }

    output {
        File de_novo_results = hailDenovoFiltering.de_novo_results
        File de_novo_vep = hailDenovoFiltering.de_novo_vep
        File de_novo_ht = hailDenovoFiltering.de_novo_ht
        File tdt_mt = hailDenovoFiltering.tdt_mt
        File tdt_parent_aware_mt = hailDenovoFiltering.tdt_parent_aware_mt
    }
}

task hailDenovoFiltering {
    input {
        File corrected_ped
        String filtered_mt
        String cohort_prefix
        String bucket_id
        String hail_denovo_filtering_script
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
        curl ~{hail_denovo_filtering_script} > hail_denovo_filtering_script.py
        python3 hail_denovo_filtering_script.py ~{filtered_mt} ~{cohort_prefix} ~{corrected_ped} ~{bucket_id} \
        ~{cpu_cores} ~{memory}
    }

    output {
        File de_novo_results = "~{cohort_prefix}_wes_final_denovo.txt"
        File de_novo_vep = "~{cohort_prefix}_wes_final_denovo_vep.txt"
        String de_novo_ht = "~{cohort_prefix}_wes_final_denovo.ht"
        String tdt_mt = "~{bucket_id}/hail/~{cohort_prefix}_trio_tdt.mt"
        String tdt_parent_aware_mt = "~{bucket_id}/hail/~{cohort_prefix}_parent_aware_trio_tdt.mt"
    }
}
