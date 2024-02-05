version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow annotateStep3 {
    input {
        Array[Array[File]]? vep_annotated_final_vcf
        Array[File]? vep_vcf_files
        Array[File] split_trio_vcfs
        File ped_uri
        String mpc_dir
        File mpc_chr22_file
        File loeuf_file
        String annotate_step03_script
        String vep_hail_docker
    }

    if (defined(vep_annotated_final_vcf)) {
        Array[File] vep_annotated_final_vcf_arr = flatten(select_first([vep_annotated_final_vcf]))
    }
    File vep_uri = select_first([vep_vcf_files, vep_annotated_final_vcf_arr])[0]

    scatter (vcf_uri in split_trio_vcfs) {
        call annotateStep03 {
            input:
                vcf_uri=vcf_uri,
                vep_uri=vep_uri,
                ped_uri=ped_uri,
                mpc_dir=mpc_dir,
                mpc_chr22_file=mpc_chr22_file,
                loeuf_file=loeuf_file,
                annotate_step03_script=annotate_step03_script,
                vep_hail_docker=vep_hail_docker
        }
    }

    output {
        Array[File] split_trio_annot_tsv = annotateStep03.split_trio_annot_tsv
    }
}

task annotateStep03 {
    input {
        File vcf_uri
        File vep_uri
        File ped_uri
        String mpc_dir
        File mpc_chr22_file
        File loeuf_file
        String annotate_step03_script
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
        curl ~{annotate_step03_script} > annotate.py
        python3.9 annotate.py ~{vcf_uri} ~{vep_uri} ~{ped_uri} ~{mpc_dir} ~{mpc_chr22_file} ~{loeuf_file} ~{cpu_cores} ~{memory}
    }

    output {
        File split_trio_annot_tsv = basename(vcf_uri, '.vcf') + "_annot.tsv"
    }
}