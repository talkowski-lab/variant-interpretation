version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow step1 {
    input {
        File vcf_file
        File corrected_ped
        File purcell5k
        File mpc_chr22_file
        String mpc_dir
        String bucket_id
        String gnomad_ht_uri
        String cohort_prefix
        String hail_annotation_script
        String vep_hail_docker
    }

    call hailAnnotate {
        input:
            vcf_file=vcf_file,
            corrected_ped=corrected_ped,
            purcell5k=purcell5k,
            mpc_chr22_file=mpc_chr22_file,
            mpc_dir=mpc_dir,
            gnomad_ht_uri=gnomad_ht_uri,
            bucket_id=bucket_id,
            cohort_prefix=cohort_prefix,
            hail_annotation_script=hail_annotation_script,
            vep_hail_docker=vep_hail_docker
    }

    output {
        File annot_mt = hailAnnotate.annot_mt
        File sample_qc_info = hailAnnotate.sample_qc_info
        File pca_score_table_5k = hailAnnotate.pca_score_table_5k
        File pca_loading_table_5k = hailAnnotate.pca_loading_table_5k
    }
}

task hailAnnotate {
    input {
        File vcf_file
        File corrected_ped
        File purcell5k
        File mpc_chr22_file
        String mpc_dir
        String bucket_id
        String gnomad_ht_uri
        String cohort_prefix
        String hail_annotation_script
        String vep_hail_docker
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

    command {
        curl ~{hail_annotation_script} > hail_annotation_script.py
        python3.9 hail_annotation_script.py ~{vcf_file} ~{cohort_prefix} ~{bucket_id} ~{corrected_ped} \
        ~{gnomad_ht_uri} ~{mpc_dir} ~{mpc_chr22_file} ~{purcell5k} ~{cpu_cores} ~{memory}
    }

    output {
        String annot_mt = "~{bucket_id}/hail/~{cohort_prefix}_wes_denovo_annot.mt"
        File sample_qc_info = "~{cohort_prefix}_wes_post_annot_sample_QC_info.txt"
        String pca_score_table_5k = "~{bucket_id}/hail/~{cohort_prefix}_wes_pca_score_table_5k.ht"
        String pca_loading_table_5k = "~{bucket_id}/hail/~{cohort_prefix}_wes_pca_loading_table_5k.ht"
    }
}
