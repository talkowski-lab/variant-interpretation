version 1.0 

workflow step5 {
    input {
        File ped_uri
        Array[Array[File]] split_trio_vcfs  # for input directory (from step 4)
        Array[Array[File]] trio_denovo_vcf  # for output directory
        File merge_vcf_to_tsv_fullQC_py
        String trio_denovo_docker
        String cohort_prefix

    }

    call merge_vcf_to_tsv_fullQC {
        input:
            ped_uri=ped_uri,
            split_trio_vcfs=split_trio_vcfs,
            trio_denovo_vcf=trio_denovo_vcf,
            merge_vcf_to_tsv_fullQC_py=merge_vcf_to_tsv_fullQC_py,
            trio_denovo_docker=trio_denovo_docker,
            cohort_prefix=cohort_prefix
    }

    output {
        File vcf_metrics_tsv = merge_vcf_to_tsv_fullQC.output_tsv
    }
}

task merge_vcf_to_tsv_fullQC {
    input {
        File ped_uri
        File merge_vcf_to_tsv_fullQC_py
        Array[Array[File]] split_trio_vcfs 
        Array[Array[File]] trio_denovo_vcf
        String trio_denovo_docker
        String cohort_prefix
    }

    runtime {
        docker: trio_denovo_docker
    }

    command {
        input_dir=$(dirname ~{split_trio_vcfs[0][0]})
        output_dir=$(dirname ~{trio_denovo_vcf[0][0]})
        python3 ~{merge_vcf_to_tsv_fullQC_py} -d $output_dir -i $input_dir -p ~{ped_uri} -o ~{cohort_prefix}_dnm.tsv
    }

    output {
        File output_tsv = cohort_prefix + '_dnm.tsv'
    }
}