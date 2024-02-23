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
        File ped_uri
        File purcell5k
        File mpc_chr22_file
        String mpc_dir
        String gnomad_ht_uri
        String cohort_prefix
        String hail_annotation_script
        String hail_docker
        String bucket_id=false
        RuntimeAttr? runtime_attr_override
    }

    if (bucket_id=='false') {
        call hailAnnotate {
            input:
                vcf_file=vcf_file,
                ped_uri=ped_uri,
                purcell5k=purcell5k,
                mpc_chr22_file=mpc_chr22_file,
                mpc_dir=mpc_dir,
                gnomad_ht_uri=gnomad_ht_uri,
                cohort_prefix=cohort_prefix,
                hail_annotation_script=hail_annotation_script,
                hail_docker=hail_docker,
                runtime_attr_override=runtime_attr_override
        }
    }

    if (bucket_id!='false') {
        call hailAnnotateRemote {
            input:
                vcf_file=vcf_file,
                ped_uri=ped_uri,
                purcell5k=purcell5k,
                mpc_chr22_file=mpc_chr22_file,
                mpc_dir=mpc_dir,
                gnomad_ht_uri=gnomad_ht_uri,
                bucket_id=bucket_id,
                cohort_prefix=cohort_prefix,
                hail_annotation_script=hail_annotation_script,
                hail_docker=hail_docker,
                runtime_attr_override=runtime_attr_override
        }
    }

    output {
        String annot_mt = select_first([hailAnnotate.annot_mt, hailAnnotateRemote.annot_mt])
        File sample_qc_info = select_first([hailAnnotate.sample_qc_info, hailAnnotateRemote.sample_qc_info])
        String pca_score_table_5k = select_first([hailAnnotate.pca_score_table_5k, hailAnnotateRemote.pca_score_table_5k])
        String pca_loading_table_5k = select_first([hailAnnotate.pca_loading_table_5k, hailAnnotateRemote.pca_loading_table_5k])
    }
}

task hailAnnotate {
    input {
        File vcf_file
        File ped_uri
        File purcell5k
        File mpc_chr22_file
        String mpc_dir
        String gnomad_ht_uri
        String cohort_prefix
        String hail_annotation_script
        String hail_docker
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

    command {
        curl ~{hail_annotation_script} > hail_annotation_script.py
        python3 hail_annotation_script.py ~{vcf_file} ~{cohort_prefix} ~{ped_uri} \
        ~{gnomad_ht_uri} ~{mpc_dir} ~{mpc_chr22_file} ~{purcell5k} ~{cpu_cores} ~{memory}
        
        tar -zcvf ~{cohort_prefix}_wes_denovo_annot.mt.gz ~{cohort_prefix}_wes_denovo_annot.mt
        tar -zcvf ~{cohort_prefix}_wes_pca_score_table_5k.ht.gz ~{cohort_prefix}_wes_pca_score_table_5k.ht
        tar -zcvf ~{cohort_prefix}_wes_pca_loading_table_5k.ht.gz ~{cohort_prefix}_wes_pca_loading_table_5k.ht
    }

    output {
        File annot_mt = "~{cohort_prefix}_wes_denovo_annot.mt.gz"
        File sample_qc_info = "~{cohort_prefix}_wes_post_annot_sample_QC_info.txt"
        File pca_score_table_5k = "~{cohort_prefix}_wes_pca_score_table_5k.ht.gz"
        File pca_loading_table_5k = "~{cohort_prefix}_wes_pca_loading_table_5k.ht.gz"
    }
}

task hailAnnotateRemote {
    input {
        File vcf_file
        File ped_uri
        File purcell5k
        File mpc_chr22_file
        String bucket_id
        String mpc_dir
        String gnomad_ht_uri
        String cohort_prefix
        String hail_annotation_script
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_file, 'GB')
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
        curl ~{hail_annotation_script} > hail_annotation_script.py
        python3 hail_annotation_script.py ~{vcf_file} ~{cohort_prefix} ~{ped_uri} \
        ~{gnomad_ht_uri} ~{mpc_dir} ~{mpc_chr22_file} ~{purcell5k} ~{cpu_cores} ~{memory} ~{bucket_id}
    }

    output {
        String annot_mt = read_lines('vcf_file.txt')[0]
        File sample_qc_info = "~{cohort_prefix}_wes_post_annot_sample_QC_info.txt"
        String pca_score_table_5k = read_lines('vcf_file.txt')[1]
        String pca_loading_table_5k = read_lines('vcf_file.txt')[2]
    }
}
