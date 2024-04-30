
version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow step4 {
    input {
        File de_novo_merged
        String cohort_prefix
        String final_filtering_script
        String hail_docker
        Int vqslod_cutoff_snv=-20
        Int vqslod_cutoff_indel=-2
        Float af_threshold=0.005
    }

    call finalFiltering {
        input:
        de_novo_merged=de_novo_merged,
        cohort_prefix=cohort_prefix,
        final_filtering_script=final_filtering_script,
        hail_docker=hail_docker,
        vqslod_cutoff_snv=vqslod_cutoff_snv,
        vqslod_cutoff_indel=vqslod_cutoff_indel,
        af_threshold=af_threshold
    }

    output {
        File de_novo_final = finalFiltering.de_novo_final
    }
}

task finalFiltering {
    input {
        File de_novo_merged
        String cohort_prefix
        String final_filtering_script
        String hail_docker
        Int vqslod_cutoff_snv=-20
        Int vqslod_cutoff_indel=-2
        Float af_threshold=0.001
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(de_novo_merged, 'GB')
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
    curl ~{final_filtering_script} > final_filtering.py
    python3 final_filtering.py ~{de_novo_merged} ~{cohort_prefix} ~{vqslod_cutoff_snv} ~{vqslod_cutoff_indel} \
        ~{af_threshold} ~{cpu_cores} ~{memory} > stdout
    >>>

    output {
        File de_novo_final = cohort_prefix + '_de_novo_filtered_final.tsv'
    }
}
