version 1.0

import "wes-denovo-helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow getBAFinBED {
    input {
        File bed_file
        Array[File] vep_vcf_files
        String cohort_prefix
        String get_baf_script
        String plot_baf_script
        String hail_docker
    }

    scatter (vep_file in vep_vcf_files) {
        call getBAF {
            input:
            bed_file=bed_file,
            vep_file=vep_file,
            cohort_prefix=cohort_prefix,
            get_baf_script=get_baf_script,
            hail_docker=hail_docker
        }
    }

    call helpers.mergeResults as mergeResults {
        input:
        tsvs=getBAF.baf_tsv,
        hail_docker=hail_docker,
        merged_filename=cohort_prefix + '_AB_SNVs_locus_intervals'
    }

    call plotBAF {
        input:
        merged_baf_tsv=mergeResults.merged_tsv,
        hail_docker=hail_docker,
        plot_baf_script=plot_baf_script
    }
    output {
        File merged_baf_tsv = mergeResults.merged_tsv
        Array[File] baf_plots = plotBAF.baf_plots
    }
}

task getBAF {
    input {
        File bed_file
        File vep_file
        String cohort_prefix
        String get_baf_script
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vep_file, "GB")
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

    String output_name = basename(vep_file, '.vcf.bgz') + '_AB_SNVs_locus_intervals.tsv'

    command <<<
    curl ~{get_baf_script} > get_baf.py
    python3 get_baf.py ~{bed_file} ~{cohort_prefix} ~{vep_file} ~{cpu_cores} ~{memory} ~{output_name} > stdout
    >>>

    output {
        File baf_tsv = output_name
    }
}

task plotBAF {
    input {
        File merged_baf_tsv
        String hail_docker
        String plot_baf_script
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(merged_baf_tsv, "GB")
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
    curl ~{plot_baf_script} > plot_baf.py
    python3 plot_baf.py ~{merged_baf_tsv}
    >>>

    output {
        Array[File] baf_plots = glob('*.png')
    }
}