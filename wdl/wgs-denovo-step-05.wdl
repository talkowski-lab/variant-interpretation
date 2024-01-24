version 1.0 

workflow step5 {
    input {
        File ped_uri
        Array[File] split_trio_vcfs  # for input directory (from step 4)
        Array[File] trio_denovo_vcf  # for output directory
        String merge_vcf_to_tsv_fullQC_script
        String trio_denovo_docker
        String cohort_prefix

    }

    call merge_vcf_to_tsv_fullQC {
        input:
            ped_uri=ped_uri,
            split_trio_vcfs=split_trio_vcfs,
            trio_denovo_vcf=trio_denovo_vcf,
            merge_vcf_to_tsv_fullQC_script=merge_vcf_to_tsv_fullQC_script,
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
        String merge_vcf_to_tsv_fullQC_script
        Array[File] split_trio_vcfs 
        Array[File] trio_denovo_vcf
        String trio_denovo_docker
        String cohort_prefix
    }

    runtime {
        docker: trio_denovo_docker
    }

    command {
        input_dir=$(dirname ~{split_trio_vcfs[0]})
        output_dir=$(dirname ~{trio_denovo_vcf[0]})
        curl ~{merge_vcf_to_tsv_fullQC_script} > merge_vcf_to_tsv_fullQC.py
        python3 merge_vcf_to_tsv_fullQC.py -d $output_dir -i $input_dir -p ~{ped_uri} -o ~{cohort_prefix}_dnm.tsv
    }

    output {
        File output_tsv = cohort_prefix + '_dnm.tsv'
    }
}