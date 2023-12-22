version 1.0 

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow filterRareVariants {
    input {
        File trio_uri
        File ped_uri
        Array[Array[File]] vep_annotated_final_vcf
        # Array[Array[File]] vep_annotated_final_vcf_idx
        String cohort_prefix
        String sv_base_mini_docker
        Float AF_threshold=0.005
        Int AC_threshold=1
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_split_vcf
        RuntimeAttr? runtime_attr_filter_vcf
    }

    scatter (cohort_vcf_files in vep_annotated_final_vcf) {
        call mergeVCFs as mergeSharded {
            input:
                vcf_files=cohort_vcf_files,
                sv_base_mini_docker=sv_base_mini_docker,
                cohort_prefix=cohort_prefix,
                merge_or_concat='concat',
                runtime_attr_override=runtime_attr_merge_vcfs
        }
    }
    call mergeVCFs as mergeCohort {
    input:
        vcf_files=mergeSharded.merged_vcf_file,
        sv_base_mini_docker=sv_base_mini_docker,
        cohort_prefix=cohort_prefix,
        merge_or_concat='concat',
        runtime_attr_override=runtime_attr_merge_vcfs
    }

    call splitTrioVCFs {
        input:
            trio_uri=trio_uri,
            vcf_file=mergeCohort.merged_vcf_file,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix,
            runtime_attr_override=runtime_attr_split_vcf
    }

    scatter (vcf_file in splitTrioVCFs.split_trio_vcfs) {
        call filterRareVariants {
            input:
                vcf_file=vcf_file,
                AC_threshold=AC_threshold,
                AF_threshold=AF_threshold,
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_filter_vcf
        }
    }

    call mergeVCFs as mergeFiltered {
        input:
            vcf_files=filterRareVariants.filtered_vcf_file,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix,
            merge_or_concat='merge',
            runtime_attr_override=runtime_attr_merge_vcfs
    }

    output {
        File merged_filtered_vcf_file=mergeFiltered.merged_vcf_file
        File merged_filtered_vcf_idx=mergeFiltered.merged_vcf_idx
    }
}

task mergeVCFs {
    input {
        Array[File] vcf_files
        String sv_base_mini_docker
        String cohort_prefix
        String merge_or_concat    
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_files, "GB")
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

    String merged_vcf_name="~{cohort_prefix}.merged.vcf.gz"

    String merge_or_concat_new = if merge_or_concat == 'concat' then 'concat -n'  else merge_or_concat

    command <<<
        set -euo pipefail
        VCFS="~{write_lines(vcf_files)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        bcftools ~{merge_or_concat_new} --no-version -Oz --file-list vcfs_sorted.list --output ~{merged_vcf_name}
        bcftools index -t ~{merged_vcf_name}
    >>>

    output {
        File merged_vcf_file=merged_vcf_name
        File merged_vcf_idx=merged_vcf_name + ".tbi"
    }
}

task splitTrioVCFs {
    input {
        File trio_uri
        File vcf_file
        String sv_base_mini_docker
        String cohort_prefix
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command {
        cat ~{trio_uri} | tail -n +2 | cut -f3-5 | tr '\t' ',' > samples.txt
        cat ~{trio_uri} | tail -n +2 | cut -f2-3 | tr '\t' '_trio_' > filenames.txt
        paste samples.txt samples.txt filenames.txt > trio.list
        bcftools +split -S trio.list -Oz -o split_trio_vcfs ~{vcf_file}
    }

    output {
        Array[File] split_trio_vcfs = glob("split_trio_vcfs/*")
    }
}

task filterRareVariants {
    input {
        File vcf_file
        Int AC_threshold
        Float AF_threshold
        String sv_base_mini_docker
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String new_filename = basename(vcf_file, '.vcf.gz') + '_ultra_rare.vcf.gz'
    
    command <<<
        bcftools view -i "(INFO/AF<~{AF_threshold} || INFO/AC=~{AC_threshold}) && GT[0]='het' && ((GT[1]!='het' && GT[2]='het') || (GT[1]='het' && GT[2]!='het'))" \
        -Oz -o ~{new_filename} ~{vcf_file}
        bcftools index -t ~{new_filename}
    >>>

    output {
        File filtered_vcf_file = new_filename
        File filtered_vcf_idx = new_filename + ".tbi"
    }
}