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
        Array[File]? vep_vcf_files
        Array[File]? vep_annotated_final_vcf
        File? merged_vep_file
        File ped_uri
        File bed_file
        File ancestry_labels_1kg
        File somalier_1kg_tar
        String correct_somalier_ped_python_script
        String cohort_prefix
        String somalier_docker
        String sv_base_mini_docker
        String hail_docker
        Boolean subset_ped=true
        Boolean infer_ped=true
        RuntimeAttr? runtime_attr_relatedness
        RuntimeAttr? runtime_attr_correct
    }

    if (!defined(merged_vep_file)) {
        Array[File] vep_files = select_first([vep_vcf_files, vep_annotated_final_vcf])
        
        scatter (vcf_uri in vep_files) {
            call subsetVCFs {
                input:
                    bed_file=bed_file,
                    vcf_uri=vcf_uri,
                    vcf_idx=vcf_uri+'.tbi',
                    somalier_docker=somalier_docker
            }
        }

        call mergeVCFs {
        input:
            vcf_files=subsetVCFs.subset_vcf,
            vcf_files_idx=subsetVCFs.subset_vcf_idx,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix,
            merge_or_concat='concat'
        }
    }

    File merged_vcf_file = select_first([merged_vep_file, mergeVCFs.merged_vcf_file])
    
    call relatedness {
        input:
            sites_uri=sites_uri,
            hg38_fasta=hg38_fasta,
            vcf_uri=merged_vcf_file,
            ped_uri=ped_uri,
            ancestry_labels_1kg=ancestry_labels_1kg,
            somalier_1kg_tar=somalier_1kg_tar,
            cohort_prefix=cohort_prefix,
            somalier_docker=somalier_docker,
            infer_ped=infer_ped,
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
        Boolean infer_ped
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

    String infer_string = if infer_ped then "--infer" else ""
    command {
        bcftools index ~{vcf_uri}
        somalier extract -d extracted/ --sites ~{sites_uri} -f ~{hg38_fasta} ~{vcf_uri}
        somalier relate ~{infer_string} --ped ~{ped_uri} -o ~{cohort_prefix} extracted/*.somalier

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
        cat ~{merged_vcf_name} | zcat | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' > ~{basename(sorted_vcf_name, '.gz')}
        bgzip ~{basename(sorted_vcf_name, '.gz')}
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
        File out_samples
        String correct_somalier_ped_python_script
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
        curl ~{correct_somalier_ped_python_script} > correct_somalier_ped_python_script.py
        python3 correct_somalier_ped_python_script.py ~{out_samples} ~{ped_uri} ~{cohort_prefix} ~{subset_ped} > stdout
    }

    output {
        File corrected_ped = cohort_prefix + "_ped_corrected.ped"
        File somalier_errors = cohort_prefix + "_somalier_errors.tsv"
    }
}