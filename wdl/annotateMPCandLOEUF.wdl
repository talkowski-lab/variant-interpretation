version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow annotateMPCandLOEUF {
    input {
        File vcf_metrics_tsv
        String mpc_dir
        File mpc_chr22_file
        File loeuf_file
        File annotate_mpc_loeuf_script
        String hail_docker
        RuntimeAttr? runtime_attr_hail_annotate
    }

    call annotateMPCandLOEUF {
        input:
            vcf_metrics_tsv=vcf_metrics_tsv,
            mpc_dir=mpc_dir,
            mpc_chr22_file=mpc_chr22_file,
            loeuf_file=loeuf_file,
            annotate_mpc_loeuf_script=annotate_mpc_loeuf_script,
            hail_docker=hail_docker,
            runtime_attr_override=runtime_attr_hail_annotate
    }

    output {
        File vcf_metrics_tsv_annot = annotateMPCandLOEUF.vcf_metrics_tsv_annot
    }
}

task annotateMPCandLOEUF {
    input {
        File vcf_metrics_tsv
        String mpc_dir
        File mpc_chr22_file
        File loeuf_file
        File annotate_mpc_loeuf_script
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([vcf_metrics_tsv, loeuf_file], "GB")
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

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    
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
        python3 ~{annotate_mpc_loeuf_script} ~{vcf_metrics_tsv} ~{mpc_dir} ~{mpc_chr22_file} ~{loeuf_file} ~{cpu_cores} ~{memory}
    }

    output {
        File vcf_metrics_tsv_annot = basename(vcf_metrics_tsv, '.tsv') + '_with_mpc_loeuf.tsv'
    }
}