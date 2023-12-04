version 1.0 

workflow step3 {
    input {
        File ped_uri
        File merged_preprocessed_vcf_file
        String hail_docker
        File uberSplit_v3_py
        Int batch_size
    }

    String cohort_prefix = basename(merged_preprocessed_vcf_file, '.vep.merged.vcf.gz')
    String stats_file = cohort_prefix + "_stats.txt"
    call uberSplit_v3 {
        input:
            ped_uri=ped_uri,
            vcf_file=merged_preprocessed_vcf_file,
            hail_docker=hail_docker,
            cohort_prefix=cohort_prefix,
            stats_file=stats_file,
            uberSplit_v3_py=uberSplit_v3_py,
            batch_size=batch_size
    }

    output {
        Array[File] split_trio_vcfs = uberSplit_v3.split_trio_vcfs
        File stats_files = uberSplit_v3.stats_file_out
    }
}

task uberSplit_v3 {
    input {
        File ped_uri
        File vcf_file
        String hail_docker
        String cohort_prefix
        String stats_file
        File uberSplit_v3_py       
        Int batch_size
    }

    runtime {
        docker: hail_docker
    }

    command {
        mkdir -p ~{cohort_prefix}
        python3 ~{uberSplit_v3_py} ~{ped_uri} ~{vcf_file} ~{cohort_prefix} ~{stats_file} ~{batch_size}
    }

    output {
        Array[File] split_trio_vcfs = glob(cohort_prefix + "/*")
        File stats_file_out = stats_file
    }
}