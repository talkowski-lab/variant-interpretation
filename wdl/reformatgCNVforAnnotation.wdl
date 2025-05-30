version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow reformatgCNV {
    input {
        File bed_file
        String prefix
        String sv_base_mini_docker
        Int num_records
        Int win_dist
        String reformat_gCNV_for_merging_script
        String bed2vcf_script
        RuntimeAttr? runtime_attr_reformat
        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_bed2vcf
        RuntimeAttr? runtime_attr_sortVcf
    }

    call reformatBed{
        input:
            bed_file=bed_file,
            sv_base_mini_docker=sv_base_mini_docker,
            prefix=prefix,
            runtime_attr_override=runtime_attr_reformat,
            reformat_gCNV_for_merging_script=reformat_gCNV_for_merging_script
    }

    call mergeBed{
        input:
            bed_file=reformatBed.reformat_bed_merge,
            sv_base_mini_docker=sv_base_mini_docker,
            prefix=prefix,
            win_dist=win_dist,
            runtime_attr_override=runtime_attr_merge
    }

    call bed2vcf{
        input:
            bed_file=reformatBed.reformat_bed,
            bed_file_merged=mergeBed.merged_bed,
            sv_base_mini_docker=sv_base_mini_docker,
            prefix=prefix,
            runtime_attr_override=runtime_attr_bed2vcf,
            bed2vcf_script=bed2vcf_script
    }

    call sortVCF{
        input:
            vcf_file=bed2vcf.vcf_file,
            sv_base_mini_docker=sv_base_mini_docker,
            prefix=prefix,
            runtime_attr_override=runtime_attr_sortVcf
    }

    output {
        File final_vcf = sortVCF.final_vcf_file
    }
}

task reformatBed {
    input {
        File bed_file
        String sv_base_mini_docker
        String prefix
        String reformat_gCNV_for_merging_script
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(bed_file, "GB")
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        curl ~{reformat_gCNV_for_merging_script} > reformat_gCNV_for_merging.R
        Rscript reformat_gCNV_for_merging.R \
            -i ~{bed_file} \
            -p ~{prefix}
    >>>

    output {
        File reformat_bed = "~{prefix}.ref.bed"
        File reformat_bed_merge = "~{prefix}.ref.for_merging.bed"
    }
}

task mergeBed {
    input {
        File bed_file
        String sv_base_mini_docker
        String prefix
        Int win_dist #default of 15000000
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(bed_file, "GB")
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        bedtools merge -d ~{win_dist} -i ~{bed_file} -c 4 -o collapse -delim ',' > ~{prefix}.merged.bed
    >>>

    output {
        File merged_bed = "~{prefix}.merged.bed"
    }
}

task bed2vcf {
    input {
        File bed_file
        File bed_file_merged
        String sv_base_mini_docker
        String prefix
        String bed2vcf_script
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(bed_file, "GB")
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail

        curl ~{bed2vcf_script} > bed2vcf.R

        Rscript bed2vcf.R \
            -i ~{bed_file} \
            -m ~{bed_file_merged} \
            -o ~{prefix}.vcf
    >>>

    output {
        File vcf_file = "~{prefix}.vcf"
    }
}

task sortVCF {
    input {
        File vcf_file
        String sv_base_mini_docker
        String prefix
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail

        vcf-sort ~{vcf_file} | bgzip -c > ~{prefix}.for_annotation.vcf.gz
        tabix -p vcf ~{prefix}.for_annotation.vcf.gz
    >>>

    output {
        File final_vcf_file = "~{prefix}.for_annotation.vcf.gz"
        File final_vcf_file_idx = "~{prefix}.for_annotation.vcf.gz.tbi"
    }
}