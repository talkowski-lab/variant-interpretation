version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow runSomalier {
    input {
        File sites_uri
        File hg38_fasta
        Array[Array[Array[File]]] vep_annotated_final_vcf
        File ped_uri
        File bed_file
        File ancestry_labels_1kg
        File correct_somalier_ped_python_script
        File somalier_1kg_tar
        String cohort_prefix
        String somalier_docker
        String sv_base_mini_docker
        String hail_docker
        Boolean subset_ped=true   
        Boolean overlapping_samples=false     
        RuntimeAttr? runtime_attr_relatedness
        RuntimeAttr? runtime_attr_correct
    }
    scatter (cohort_vcf_files in vep_annotated_final_vcf) {
        scatter (vcf_files in cohort_vcf_files) {
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
    }

    if (overlapping_samples) {
        call splitMergeVCFs {
            input:
                vcf_files=mergeCohort.merged_vcf_file,
                vcf_files_idx=mergeCohort.merged_vcf_idx,
                sv_base_mini_docker=sv_base_mini_docker,
                cohort_prefix=cohort_prefix
        }
    }
    if (!overlapping_samples) {
        call mergeVCFs {
            input:
                vcf_files=mergeCohort.merged_vcf_file,
                vcf_files_idx=mergeCohort.merged_vcf_idx,
                sv_base_mini_docker=sv_base_mini_docker,
                cohort_prefix=cohort_prefix,
                merge_or_concat='merge'
        }
    }

    call relatedness {
        input:
            sites_uri=sites_uri,
            hg38_fasta=hg38_fasta,
            vcf_uri=select_first([splitMergeVCFs.merged_vcf_file, mergeVCFs.merged_vcf_file]),
            ped_uri=ped_uri,
            ancestry_labels_1kg=ancestry_labels_1kg,
            somalier_1kg_tar=somalier_1kg_tar,
            cohort_prefix=cohort_prefix,
            somalier_docker=somalier_docker,
            runtime_attr_override=runtime_attr_relatedness
    }

    call correctPedigree {
        input:
            ped_uri=ped_uri,
            correct_somalier_ped_python_script=correct_somalier_ped_python_script,
            out_samples=relatedness.out_samples,
            cohort_prefix=cohort_prefix,
            hail_docker=hail_docker,
            subset_ped=subset_ped,
            runtime_attr_override=runtime_attr_correct
    }

    output {
        File out_samples = relatedness.out_samples
        File out_pairs = relatedness.out_pairs
        File out_groups = relatedness.out_groups
        File out_html = relatedness.out_html
        File ancestry_html = relatedness.ancestry_html
        File ancestry_out = relatedness.ancestry_out
        File corrected_ped = correctPedigree.corrected_ped
        File somalier_errors = correctPedigree.somalier_errors
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

task relatedness {
    input {
        File sites_uri
        File hg38_fasta
        File vcf_uri
        File ped_uri
        File ancestry_labels_1kg
        File somalier_1kg_tar
        String cohort_prefix
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
        bcftools index ~{vcf_uri}
        somalier extract -d extracted/ --sites ~{sites_uri} -f ~{hg38_fasta} ~{vcf_uri}
        somalier relate --infer --ped ~{ped_uri} -o ~{cohort_prefix} extracted/*.somalier

        tar -xf ~{somalier_1kg_tar}
        somalier ancestry -o ~{cohort_prefix} --labels ~{ancestry_labels_1kg} 1kg-somalier/*.somalier ++ extracted/*.somalier
    }

    output {
        File out_samples = cohort_prefix + ".samples.tsv" # creates a .ped like vep_annotated_final_vcf with extra QC columns
        File out_pairs = cohort_prefix + ".pairs.tsv" # shows IBS for all possible sample pairs
        File out_groups = cohort_prefix + ".groups.tsv" # shows pairs of samples above a certain relatedness
        File out_html = cohort_prefix + ".html" # interactive html
        File ancestry_html = cohort_prefix + ".somalier-ancestry.html"
        File ancestry_out = cohort_prefix + ".somalier-ancestry.tsv"
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

task correctPedigree {
    input {
        File ped_uri
        File correct_somalier_ped_python_script
        File out_samples
        String cohort_prefix
        String hail_docker
        Boolean subset_ped=true
        RuntimeAttr? runtime_attr_override    
    }

    Float relatedness_size = size(ped_uri, "GB") 
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
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command {
        python3 ~{correct_somalier_ped_python_script} ~{out_samples} ~{ped_uri} ~{cohort_prefix} ~{subset_ped} > stdout
    }

    output {
        File corrected_ped = cohort_prefix + "_ped_corrected.ped"
        File somalier_errors = cohort_prefix + "_somalier_errors.tsv"
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

            shared_samples=$(comm -12 merged_vcf_samples.txt vcf_samples.txt | tr "\n" "," | head -c -2)
            bcftools view -s $shared_samples -Oz -o tmp_merged_shared.vcf.gz $merged_vcf
            bcftools view -s $shared_samples -Oz -o tmp_shared.vcf.gz $cur_vcf
            bcftools concat -n --no-version -Oz --output tmp_merged.vcf.gz tmp_merged_shared.vcf.gz tmp_shared.vcf.gz 

            bcftools view -s "^$shared_samples" -Oz -o tmp_merged_unique.vcf.gz $merged_vcf
            bcftools view -s "^$shared_samples" -Oz -o tmp_unique.vcf.gz $cur_vcf
            bcftools merge --no-version -Oz --output tmp_unique.vcf.gz tmp_merged_unique.vcf.gz tmp_unique.vcf.gz 

            bcftools merge --no-version -Oz --output $new_merged_vcf tmp_merged.vcf.gz tmp_unique.vcf.gz
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
