version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow annotateStep1 {
    input {
        Array[Array[File]]? vep_annotated_final_vcf
        Array[File]? vep_vcf_files
        File merged_preprocessed_vcf_file
        String mpc_dir
        File mpc_chr22_file
        File loeuf_file
        String cohort_prefix
        String annotate_vcf_script
        String vep_hail_docker
    }

    if (defined(vep_annotated_final_vcf)) {
        Array[File] vep_annotated_final_vcf_arr = flatten(select_first([vep_annotated_final_vcf]))
    }
    File vep_uri = select_first([vep_vcf_files, vep_annotated_final_vcf_arr])[0]

    call annotateStep01 {
        input:
            vcf_uri=merged_preprocessed_vcf_file,
            vep_uri=vep_uri,
            mpc_dir=mpc_dir,
            mpc_chr22_file=mpc_chr22_file,
            loeuf_file=loeuf_file,
            file_ext='.vcf' + sub(basename(merged_preprocessed_vcf_file), '.*.vcf', ''),
            sample='',
            annotate_vcf_script=annotate_vcf_script,
            vep_hail_docker=vep_hail_docker
    }

    output {
        File merged_preprocessed_vcf_file_annot = annotateStep01.merged_preprocessed_vcf_file_annot
    }
}

task annotateStep01 {
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
        File merged_preprocessed_vcf_file_annot = basename(vcf_uri, '.vcf') + "_annot.tsv"
    }
}