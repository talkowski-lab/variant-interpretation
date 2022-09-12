version 1.0
    
workflow deNovoSV {

    input {
        String bed_dir
        File ped_input
        String vcf_dir
        File disorder_input
        Array[String] contigs
        String raw_dir
        File python_config
        String variant_interpretation_docker
    }

    scatter (contig in contigs) {

        File bed_input = "~{bed_dir}/phase2.annotated.~{contig}.bed.gz"
        File vcf_input = "~{vcf_dir}/phase2.annotated.~{contig}.noheader.vcf.gz"
        File raw_input = "~{raw_dir}/~{contig}_all_batches_m01_ref_sort.bed.gz"

        call getDeNovo{
            input:
                bed_input=bed_input,
                ped_input=ped_input,
                vcf_input=vcf_input,
                disorder_input=disorder_input,
                chromosome=contig,
                raw_input=raw_input,
                python_config=python_config,
                variant_interpretation_docker=variant_interpretation_docker
        }
    }
}

task getDeNovo{
    input{
        File bed_input
        File ped_input
        File vcf_input
        File disorder_input
        String chromosome
        File raw_input
        File python_config
        String variant_interpretation_docker
    }

    output{
        File denovo_output = "~{chromosome}.denovo.bed"
        File denovo_outliers = "~{chromosome}.denovo.outliers.bed"
    }

    command <<<
            python /src/variant-interpretation/scripts/deNovoSVs.py \
                --bed ~{bed_input} \
                --ped ~{ped_input} \
                --vcf ~{vcf_input} \
                --disorder ~{disorder_input} \
                --out ~{chromosome}.denovo.bed \
                --raw ~{raw_input} \
                --config ~{python_config} \
                --verbose True \
                --outliers ~{chromosome}.denovo.outliers.bed
    >>>

    runtime {
        memory: "12 GiB"
        disks: "local-disk 4 HDD"
        cpu: 1
        preemptible: 3
        maxRetries: 1
        docker: variant_interpretation_docker
        bootDiskSizeGb: 8
    }
}
