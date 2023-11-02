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
                hail_docker=hail_docker
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
    }

    runtime {
        docker: hail_docker
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