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
        Array[Array[Array[File]]] vep_annotated_final_vcf
        File purcell5k
        String relatedness_docker
        String somalier_docker
        String sv_base_mini_docker
        String cohort_prefix
        RuntimeAttr? runtime_attr_relatedness  
        RuntimeAttr? runtime_attr_sex  
    }
    scatter (cohort_vcf_files in vep_annotated_final_vcf) {
        scatter (vcf_files in cohort_vcf_files) {
            scatter (vcf_uri in vcf_files) {
                call subsetVCFs {
                    input:
                        bed_file=purcell5k,
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
    }

    call splitMergeVCFs {
        input:
            vcf_files=mergeCohort.merged_vcf_file,
            vcf_files_idx=mergeCohort.merged_vcf_idx,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix
    }

    call VCFToolsRelatedness {
        input:
            vcf_file=splitMergeVCFs.merged_vcf_file,
            purcell5k=purcell5k,
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
        bcftools view -R ~{bed_file} ~{vcf_uri} -o ~{basename(vcf_uri, '.vcf.gz')+".somalier.subset.vcf.gz"}
        bcftools index -t ~{basename(vcf_uri, '.vcf.gz')+".somalier.subset.vcf.gz"}
    }

    output {
        File subset_vcf = basename(vcf_uri, '.vcf.gz')+".somalier.subset.vcf.gz"
        File subset_vcf_idx = basename(vcf_uri, '.vcf.gz')+".somalier.subset.vcf.gz.tbi"
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
    Float base_mem_gb = 2.0
    Float input_mem_scale = 3.0
    Float input_disk_scale = 5.0
    
    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb + input_size * input_mem_scale,
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

task splitMergeVCFs {
    input {
        Array[File] vcf_files
        Array[File] vcf_files_idx
        String sv_base_mini_docker
        String cohort_prefix
        RuntimeAttr? runtime_attr_override
    }

    #  generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
    #  generally assume working memory is ~3 * inputs
    #  CleanVcf5.FindRedundantMultiallelics
    Float input_size = size(vcf_files, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 2.0
    Float input_mem_scale = 3.0
    Float input_disk_scale = 5.0
    
    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb + input_size * input_mem_scale,
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

    command <<<
        set -euo pipefail
        VCFS="~{write_lines(vcf_files)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        vcfs=( $(cat vcfs_sorted.list) )

        for i in "${!vcfs[@]}"; 
        do     
            if [[ "$i" -eq "0" ]]; then  # first in list
                merged_vcf=${vcfs[$i]}
                continue
            fi 
            if [[ "$i" -eq "${#vcfs[@]}-1" ]]; then  # last
                new_merged_vcf=~{merged_vcf_name}
            else 
                new_merged_vcf=merged_$i.vcf.gz
            fi 
            cur_vcf=${vcfs[$i]}
            bcftools query -l $merged_vcf | sort -u > merged_vcf_samples.txt
            bcftools query -l $cur_vcf | sort -u > vcf_samples.txt

            shared_samples=$(comm -12 merged_vcf_samples.txt vcf_samples.txt | tr "\n" "," | head -c -1)
            if [ -z "$shared_samples" ]; then
                bcftools view -G -Oz -o tmp_merged_shared.vcf.gz $merged_vcf
                bcftools view -G -Oz -o tmp_shared.vcf.gz $cur_vcf
            else
                bcftools view -s $shared_samples -Oz -o tmp_merged_shared.vcf.gz $merged_vcf
                bcftools view -s $shared_samples -Oz -o tmp_shared.vcf.gz $cur_vcf
            fi
            bcftools index -t tmp_merged_shared.vcf.gz
            bcftools index -t tmp_shared.vcf.gz
            bcftools concat -a --no-version -Oz --output tmp_merged.vcf.gz tmp_merged_shared.vcf.gz tmp_shared.vcf.gz 
            bcftools sort tmp_merged.vcf.gz -o tmp_merged_sorted.vcf.gz
            bcftools index -t tmp_merged_sorted.vcf.gz

            merged_vcf_samples=$(comm -23 merged_vcf_samples.txt vcf_samples.txt | tr "\n" "," | head -c -1)
            if [ -z "$merged_vcf_samples" ]; then
                bcftools view -G -Oz -o tmp_merged_unique.vcf.gz $merged_vcf
            else
                bcftools view -s $merged_vcf_samples -Oz -o tmp_merged_unique.vcf.gz $merged_vcf
            fi
            bcftools index -t tmp_merged_unique.vcf.gz
            cur_vcf_samples=$(comm -13 merged_vcf_samples.txt vcf_samples.txt | tr "\n" "," | head -c -1)
            if [ -z "$cur_vcf_samples" ]; then
                bcftools view -G -Oz -o tmp_unique.vcf.gz $cur_vcf
            else
                bcftools view -s $cur_vcf_samples -Oz -o tmp_unique.vcf.gz $cur_vcf
            fi
            bcftools index -t tmp_unique.vcf.gz
            bcftools merge --no-version -Oz --output tmp_merged2.vcf.gz tmp_merged_unique.vcf.gz tmp_unique.vcf.gz 
            bcftools sort tmp_merged2.vcf.gz -o tmp_merged2_sorted.vcf.gz
            bcftools index -t tmp_merged2_sorted.vcf.gz
            bcftools merge --no-version -Oz --output merged.vcf.gz tmp_merged_sorted.vcf.gz tmp_merged2_sorted.vcf.gz
            bcftools sort merged.vcf.gz -o $new_merged_vcf
            bcftools index -t $new_merged_vcf
            merged_vcf=$new_merged_vcf
        done

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
        File purcell5k
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
        vcftools --gzvcf ~{vcf_file} --relatedness2 --bed ~{purcell5k}
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