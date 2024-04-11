version 1.0 

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow step3 {
    input {
        File trio_uri
        File ped_uri
        File merged_preprocessed_vcf_file
        String hail_docker
        String cohort_prefix
        String trio_denovo_docker
        String uberSplit_v3_script
        String subset_ped_script
        Int batch_size
        Boolean subset_ped=true
    }

    String stats_file = cohort_prefix + "_stats.txt"

    if (subset_ped) {
        call subsetPed {
            input:
                ped_uri=ped_uri,
                vcf_file=merged_preprocessed_vcf_file,
                trio_denovo_docker=trio_denovo_docker,
                subset_ped_script=subset_ped_script
        }
    }

    File new_ped_uri = select_first([subsetPed.new_ped_uri, ped_uri])

    call uberSplit_v3 {
        input:
            ped_uri=new_ped_uri,
            vcf_file=merged_preprocessed_vcf_file,
            hail_docker=hail_docker,
            cohort_prefix=cohort_prefix,
            stats_file=stats_file,
            uberSplit_v3_script=uberSplit_v3_script,
            batch_size=batch_size
    }

    output {
        Array[File] split_trio_vcfs = uberSplit_v3.split_trio_vcfs
        File stats_files = uberSplit_v3.stats_file_out
    }
}

task subsetPed {
    input {
        File ped_uri
        File vcf_file
        String subset_ped_script
        String trio_denovo_docker
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
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: trio_denovo_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        bcftools query -l ~{vcf_file} > samples.txt
        curl ~{subset_ped_script} > subset_ped_script.py
        python3 subset_ped_script.py samples.txt ~{ped_uri} > stdout
    >>>

    output {
        File new_ped_uri = basename(ped_uri, '.ped')+'_subset.ped'
    }
}

task splitTrioVCFs {
    input {
        File trio_uri
        File vcf_file
        File vep_annotated_final_vcf_single
        String sv_base_mini_docker
        String cohort_prefix
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB") + size(vep_annotated_final_vcf_single, "GB")
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
        bcftools head ~{vep_annotated_final_vcf_single} > og_header.txt
        grep "FILTER=" og_header.txt > new_header.txt
        bcftools annotate -h new_header.txt -Oz -o temp.vcf.gz ~{vcf_file}

        cat ~{trio_uri} | tail -n +2 | cut -f3-5 | tr '\t' ',' > samples.txt
        cat ~{trio_uri} | tail -n +2 | cut -f2-3 | tr '\t' '_trio_' > filenames.txt
        paste samples.txt samples.txt filenames.txt > trio.list
        bcftools +split -S trio.list -Ov -o split_trio_vcfs temp.vcf.gz
    }

    output {
        Array[File] split_trio_vcfs = glob("split_trio_vcfs/*")
    }
}

task uberSplit_v3 {
    input {
        File ped_uri
        File vcf_file
        String hail_docker
        String cohort_prefix
        String stats_file
        String uberSplit_v3_script       
        Int batch_size
    }

    runtime {
        docker: hail_docker
    }

    command {
        set -eou pipefail
        mkdir -p ~{cohort_prefix}
        curl ~{uberSplit_v3_script} > uberSplit_v3.py
        python3 uberSplit_v3.py ~{ped_uri} ~{vcf_file} ~{cohort_prefix} ~{stats_file} ~{batch_size}
    }

    output {
        Array[File] split_trio_vcfs = glob(cohort_prefix + "/*")
        File stats_file_out = stats_file
    }
}