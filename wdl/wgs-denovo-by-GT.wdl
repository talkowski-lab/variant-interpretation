version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow getDenovoByGT {
    input {
        File vcf_file
        File ped_uri
        Float af_threshold=0.01
        String cohort_prefix
        String vep_hail_docker
        String denovo_snv_indels_gt_script
    }

    call denovoByGT {
        input:
            vcf_file=vcf_file,
            ped_uri=ped_uri,
            af_threshold=af_threshold,
            cohort_prefix=cohort_prefix,
            vep_hail_docker=vep_hail_docker,
            denovo_snv_indels_gt_script=denovo_snv_indels_gt_script
    }

    output {
        File denovo_gt = denovoByGT.denovo_gt
    }
}

task denovoByGT {
    input {
        File vcf_file
        File ped_uri
        Float af_threshold
        String cohort_prefix
        String vep_hail_docker
        String denovo_snv_indels_gt_script
        RuntimeAttr? runtime_attr_override
    }

        Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 10.0
    RuntimeAttr runtime_default = object {
        mem_gb: 8,
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
        curl ~{denovo_snv_indels_gt_script} > denovo_snv_indels.py
        python3.9 denovo_snv_indels.py ~{vcf_file} ~{ped_uri} ~{af_threshold} ~{cohort_prefix} ~{cpu_cores} ~{memory}
    }

    output {
        File denovo_gt = "~{cohort_prefix}_denovo_GT_AF_filter.tsv.gz"
    }
}