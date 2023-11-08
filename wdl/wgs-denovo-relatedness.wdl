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
        # file can be a list of vcf files or just one vcf file
        File file
        File python_relatedness_script
        File lcr_uri
        File ped_uri
        File purcell5k
        String bucket_id
        String cohort_prefix
        String hail_docker    
        String sv_base_mini_docker    
        RuntimeAttr? runtime_attr_relatedness
    }

    String filename = basename(file)
    # if file is vcf.gz (just one file)
    Array[File] vcf_files = if (sub(filename, ".vcf.gz", "") != filename) then [file] else read_lines(file)

    File meta_uri = "~{bucket_id}/resources/metadata/~{cohort_prefix}_sample_list.txt"
    File trio_uri = "~{bucket_id}/resources/metadata/~{cohort_prefix}_trio_list.txt"
    
    scatter (vcf_uri in vcf_files) {
        call subsetVCF {
            input:
                vcf_uri=vcf_uri,
                purcell5k=purcell5k,
                cohort_prefix=cohort_prefix,
                sv_base_mini_docker=sv_base_mini_docker
        }

        call runRelatedness {
            input:
                python_relatedness_script=python_relatedness_script,
                vcf_uri=vcf_uri,
                purcell5k_subset_vcf=subsetVCF.purcell5k_subset_vcf,
                lcr_uri=lcr_uri,
                ped_uri=ped_uri,
                meta_uri=meta_uri,
                trio_uri=trio_uri,
                cohort_prefix=cohort_prefix,
                hail_docker=hail_docker,
                runtime_attr_override=runtime_attr_relatedness
        }
    }

    output {
        # Array[File] imputed_sex_res = runRelatedness.imputed_sex_res
        # Array[File] pc_relate_res = runRelatedness.pc_relate_res
        Array[File] king_res = runRelatedness.king_res
    }
}

task subsetVCF {
    input {
        File vcf_uri
        File purcell5k
        String cohort_prefix
        String sv_base_mini_docker 
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command {
        bcftools index ~{vcf_uri}
        bcftools view -R ~{purcell5k} ~{vcf_uri} -Oz -o ~{cohort_prefix}_purcell5k.vcf.gz
    }

    output {
        File purcell5k_subset_vcf = cohort_prefix + "_purcell5k.vcf.gz"
    }
}

task runRelatedness {
    input {
        File python_relatedness_script
        File purcell5k_subset_vcf
        File vcf_uri
        File lcr_uri
        File ped_uri
        File meta_uri
        File trio_uri
        String cohort_prefix
        String hail_docker 
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
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command {
        python3 ~{python_relatedness_script} ~{lcr_uri} ~{ped_uri} ~{meta_uri} ~{trio_uri} ~{purcell5k_subset_vcf}
    }

    output {
        # File imputed_sex_res = basename(vcf_uri, '.vcf.gz') + '_imputed_sex_res.tsv'
        # File pc_relate_res = basename(vcf_uri, '.vcf.gz') + '_relatedness_pc_relate_res.tsv'
        File king_res = basename(purcell5k_subset_vcf, '.vcf.gz') + '_relatedness_king_res.tsv'
    }
}