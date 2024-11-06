version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

task annotateVCF {
    input {
        File vcf_uri
        File vep_uri
        File header_file
        String mpc_dir
        File mpc_chr22_file
        File loeuf_file
        String file_ext
        String sample
        String annotate_vcf_script
        String hail_docker
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
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        curl ~{annotate_vcf_script} > annotate_vcf.py
        python3 annotate_vcf.py ~{vcf_uri} ~{vep_uri} ~{mpc_dir} ~{mpc_chr22_file} ~{loeuf_file} ~{header_file} \
        ~{file_ext} ~{sample} ~{cpu_cores} ~{memory}
    >>>

    output {
        File annotated_tsv = basename(vcf_uri, file_ext) + "_annot.tsv"
    }
}