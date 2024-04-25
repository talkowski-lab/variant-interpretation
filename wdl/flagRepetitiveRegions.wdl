version 1.0
    
struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow flagRepetitiveRegions {
    input {
        File tsv
        File repetitive_regions_bed
        String tsv_to_bed_script
        String sv_base_mini_docker
        String hail_docker
    }

    call convertTSVtoBED {
        input:
            tsv=tsv,
            tsv_to_bed_script=tsv_to_bed_script,
            hail_docker=hail_docker
    }

    call intersectBED {
        input:
            repetitive_regions_bed=repetitive_regions_bed,
            input_bed=convertTSVtoBED.bed_file,
            sv_base_mini_docker=sv_base_mini_docker
    }

    output {
        File output_bed = intersectBED.output_bed
    }
}

task convertTSVtoBED {
    input {
        File tsv
        String tsv_to_bed_script
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(tsv, "GB")
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

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command {
        curl ~{tsv_to_bed_script} > tsv_to_bed.py
        python3 tsv_to_bed.py ~{tsv}
    }
    String file_ext = if sub(basename(tsv), '\.gz', '')==basename(tsv) then '.tsv' else '.tsv.gz'
    output {
        File bed_file = basename(tsv, file_ext) + '.bed'
    }
}

task intersectBED {
    input {
        File repetitive_regions_bed
        File input_bed
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([repetitive_regions_bed, input_bed], "GB")
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

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command {
        zcat ~{repetitive_regions_bed} | bedtools intersect -a ~{input_bed} -b stdin > ~{basename(input_bed, '.bed')}_repetitive_regions.bed
    }

    output {
        File output_bed = basename(input_bed, '.bed') + '_repetitive_regions.bed'
    }
}