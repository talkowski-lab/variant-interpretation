version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow relatedness {
    input {
        Array[Array[File]] vep_annotated_final_vcf
        File bed_file
        String relatedness_docker
        String somalier_docker
        String sv_base_mini_docker
        String cohort_prefix
        RuntimeAttr? runtime_attr_relatedness  
        RuntimeAttr? runtime_attr_sex  
    }
    scatter (vcf_files in vep_annotated_final_vcf) {
        scatter (vcf_uri in vcf_files) {
            call subsetVCFs {
                input:
                    bed_file=bed_file,
                    vcf_uri=vcf_uri,
                    vcf_idx=vcf_uri+'.tbi',
                    somalier_docker=somalier_docker
            }
        }

        call mergeVCFs as mergeSharded {
            input:
                vcf_files=subsetVCFs.subset_vcf,
                vcf_files_idx=subsetVCFs.subset_vcf_idx,
                sv_base_mini_docker=sv_base_mini_docker,
                cohort_prefix=cohort_prefix,
                merge_or_concat='concat'
        }
    }
    call mergeVCFs as mergeCohort {
    input:
        vcf_files=mergeSharded.merged_vcf_file,
        vcf_files_idx=mergeSharded.merged_vcf_idx,
        sv_base_mini_docker=sv_base_mini_docker,
        cohort_prefix=cohort_prefix,
        merge_or_concat='concat'
    }

    call VCFToolsRelatedness {
        input:
            vcf_file=mergeCohort.merged_vcf_file,
            relatedness_docker=relatedness_docker,
            runtime_attr_override=runtime_attr_relatedness
    }

    # call ImputeSexPLINK {
    #     input:
    #         vcf_file=vcf_file,
    #         relatedness_docker=relatedness_docker,
    #         runtime_attr_override=runtime_attr_sex
    # }

    output {
        File out_relatedness = VCFToolsRelatedness.out_relatedness
        # File out_sex = ImputeSexPLINK.out_sex
    }
}

task subsetVCFs {
    input {
        File vcf_uri
        File vcf_idx
        File bed_file
        String somalier_docker
        RuntimeAttr? runtime_attr_override
    }

    Float relatedness_size = size(vcf_uri, "GB") 
    Float base_disk_gb = 10.0
    RuntimeAttr runtime_default = object {
                                      mem_gb: 16,
                                      disk_gb: ceil(base_disk_gb + (relatedness_size) * 5.0),
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
        docker: somalier_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command {
        bcftools view -R ~{bed_file} ~{vcf_uri} -o ~{basename(vcf_uri, '.vcf.gz')+".subset.vcf.gz"}
        bcftools index -t ~{basename(vcf_uri, '.vcf.gz')+".subset.vcf.gz"}
    }

    output {
        File subset_vcf = basename(vcf_uri, '.vcf.gz')+".subset.vcf.gz"
        File subset_vcf_idx = basename(vcf_uri, '.vcf.gz')+".subset.vcf.gz.tbi"
    }
}

task mergeVCFs {
    input {
        Array[File] vcf_files
        Array[File] vcf_files_idx
        String sv_base_mini_docker
        String cohort_prefix
        String merge_or_concat    
        RuntimeAttr? runtime_attr_override
    }

    #  generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
    #  generally assume working memory is ~3 * inputs
    #  CleanVcf5.FindRedundantMultiallelics
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
    String sorted_vcf_name="~{cohort_prefix}.merged.sorted.vcf.gz"

    String merge_or_concat_new = if merge_or_concat == 'concat' then 'concat -n'  else merge_or_concat

    command <<<
        set -euo pipefail
        VCFS="~{write_lines(vcf_files)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        bcftools ~{merge_or_concat_new} --no-version -Oz --file-list vcfs_sorted.list --output ~{merged_vcf_name}
        bcftools sort ~{merged_vcf_name} --output ~{sorted_vcf_name}
        bcftools index -t ~{sorted_vcf_name}
    >>>

    output {
        File merged_vcf_file=sorted_vcf_name
        File merged_vcf_idx=sorted_vcf_name + ".tbi"
    }
}

task VCFToolsRelatedness {
    input {
        File vcf_file
        String relatedness_docker
        RuntimeAttr? runtime_attr_override  
    }

    Float relatedness_size = size(vcf_file, "GB") 
    Float base_disk_gb = 10.0
    RuntimeAttr runtime_default = object {
                                      mem_gb: 16,
                                      disk_gb: ceil(base_disk_gb + (relatedness_size) * 5.0),
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
        docker: relatedness_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    output{
        File out_relatedness = "out.relatedness2"
    }

    command <<<
        ##Run VCFtools
        vcftools --gzvcf ~{vcf_file} --relatedness2 
    >>>
}

task ImputeSexPLINK {
    input {
        File vcf_file
        String relatedness_docker
        RuntimeAttr? runtime_attr_override  
    }

    Float relatedness_size = size(vcf_file, "GB") 
    Float base_disk_gb = 10.0
    RuntimeAttr runtime_default = object {
                                      mem_gb: 16,
                                      disk_gb: ceil(base_disk_gb + (relatedness_size) * 5.0),
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
        docker: relatedness_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command {
        plink --vcf ~{vcf_file} --split-x 'hg38' --make-bed --out ~{basename(vcf_file, '.vcf.gz') + '_split_x'}
        plink --bfile ~{basename(vcf_file, '.vcf.gz') + '_split_x'} --recode vcf bgz --out ~{basename(vcf_file, '.vcf.gz') + '_split_x'}
        plink --vcf ~{basename(vcf_file, '.vcf.gz') + '_split_x.vcf.gz'} --check-sex
    }

    output {
        File out_sex = 'plink.sexcheck'
    }
}