version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow annotateStep4 {
    input {
        Array[File]? vep_annotated_final_vcf
        Array[File]? vep_vcf_files
        Array[File] split_trio_vcfs
        String mpc_dir
        File mpc_chr22_file
        File loeuf_file
        String cohort_prefix
        String annotate_vcf_script
        String vep_hail_docker
    }

    File vep_uri = select_first([vep_vcf_files, vep_annotated_final_vcf])[0]

    scatter (vcf_uri in split_trio_vcfs) {
        call annotateStep04 {
            input:
                vcf_uri=vcf_uri,
                vep_uri=vep_uri,
                mpc_dir=mpc_dir,
                mpc_chr22_file=mpc_chr22_file,
                loeuf_file=loeuf_file,
                file_ext='.denovos.vcf.gz',
                sample=sub(basename(vcf_uri, '.denovos.vcf.gz'), '.*_trio_', ''),
                annotate_vcf_script=annotate_vcf_script,
                vep_hail_docker=vep_hail_docker
        }
    }

    call mergeResults {
        input:
            split_trio_annot_tsvs=annotateStep04.split_trio_annot_tsv_,
            cohort_prefix=cohort_prefix,
            vep_hail_docker=vep_hail_docker
    }

    output {
        File split_trio_annot_tsv = mergeResults.split_trio_annot_tsv
    }
}

task annotateStep04 {
    input {
        File vcf_uri
        File vep_uri
        String mpc_dir
        File mpc_chr22_file
        File loeuf_file
        String file_ext
        String sample
        String annotate_vcf_script
        String vep_hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_uri, "GB")
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
        curl ~{annotate_vcf_script} > annotate_vcf.py
        python3.9 annotate_vcf.py ~{vcf_uri} ~{vep_uri} ~{mpc_dir} ~{mpc_chr22_file} ~{loeuf_file} \
        ~{file_ext} ~{sample} ~{cpu_cores} ~{memory}
    }

    output {
        File split_trio_annot_tsv_ = basename(vcf_uri, '.vcf') + "_annot.tsv"
    }
}

task mergeResults {
    input {
        Array[File] split_trio_annot_tsvs
        String vep_hail_docker
        String cohort_prefix
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(split_trio_annot_tsvs, "GB")
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

    command <<<
        head -n 1 ~{split_trio_annot_tsvs[0]} > "~{cohort_prefix}_split_trio_vcfs_annot.tsv"; 
        tail -n +2 -q ~{sep=' ' split_trio_annot_tsvs} >> "~{cohort_prefix}_split_trio_vcfs_annot.tsv"
    >>>

    output {
        File split_trio_annot_tsv = cohort_prefix + '_split_trio_vcfs_annot.tsv'
    }
}