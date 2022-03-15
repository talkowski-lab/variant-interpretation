version 1.0
    
workflow deNovoSV {

    input {
        File bed_input
        File ped_input
        File vcf_input
        File disorder_input
        Array[String] contigs
        File raw_input
        String variant_interpretation_docker
    }

    scatter (contig in contigs) {
        call getDeNovo{
            input:
                bed_input=bed_input,
                ped_input=ped_input,
                vcf_input=vcf_input,
                disorder_input=disorder_input,
                chromosome=contig,
                raw_input=raw_input,
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
        String variant_interpretation_docker
    }

    output{
        File denovo_output = "~{chromosome}.denovo.bed"
        File denovo_outliers = "~{chromosome}.denovo.outliers.bed"
    }

    command <<<
            python denovo.py \
                --bed ~{bed_input} \
                --ped ~{ped_input} \
                --vcf ~{vcf_input} \
                --disorder ~{disorder_input} \
                --out ~{chromosome}.denovo.bed \
                --raw ~{raw_input} \
                --verbose True \
                --outliers ~{chromosome}.denovo.outliers.bed
    >>>

    runtime {
        memory: "24 GiB"
        disks: "local-disk 32 HDD"
        cpu: 1
        preemptible: 3
        maxRetries: 1
        docker: variant_interpretation_docker
        bootDiskSizeGb: 32
    }
}