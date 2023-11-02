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
        String bucket_id
        String cohort_prefix
        String hail_docker        
        RuntimeAttr? runtime_attr_relatedness
    }

    String filename = basename(file)
    # if file is vcf.gz (just one file)
    Array[File] vcf_files = if (sub(filename, ".vcf.gz", "") != filename) then [file] else read_lines(file)

    File meta_uri = "~{bucket_id}/resources/metadata/~{cohort_prefix}_sample_list.txt"
    File trio_uri = "~{bucket_id}/resources/metadata/~{cohort_prefix}_trio_list.txt"
    
    scatter (vcf_uri in vcf_files) {
        call runRelatedness {
            input:
                python_relatedness_script=python_relatedness_script,
                vcf_uri=vcf_uri,
                lcr_uri=lcr_uri,
                ped_uri=ped_uri,
                meta_uri=meta_uri,
                trio_uri=trio_uri,
                cohort_prefix=cohort_prefix,
                hail_docker=hail_docker,
                runtime_attr_override=runtime_attr_normalize
        }
    }

    output {
        Array[File] imputed_sex_res = runRelatedness.imputed_sex_res
        Array[File] pc_relate_res = runRelatedness.pc_relate_res
        Array[File] king_res = runRelatedness.king_res
    }
}

task runRelatedness {
    input {
        File python_relatedness_script
        File vcf_uri
        File lcr_uri
        File ped_uri
        File meta_uri
        File trio_uri
        String cohort_prefix
        String hail_docker 
        RuntimeAttr? runtime_attr_override       
    }

    RuntimeAttr runtime_default = object {
                                      mem_gb: 16,
                                      disk_gb: ceil(base_disk_gb + (vep_annotate_sizes + norm_vcf_sizes) * 5.0),
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
        python3 ~{python_relatedness_script} ~{lcr_uri} ~{ped_uri} ~{meta_uri} ~{trio_uri} ~{vcf_uri}
    }

    output {
        File imputed_sex_res = basename(vcf_uri, '.vcf.gz') + '_imputed_sex_res.tsv'
        File pc_relate_res = basename(vcf_uri, '.vcf.gz') + '_relatedness_pc_relate_res.tsv'
        File king_res = basename(vcf_uri, '.vcf.gz') + '_relatedness_king_res.tsv'
    }
}