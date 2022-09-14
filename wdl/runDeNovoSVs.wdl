version 1.0
    
import "Structs.wdl"

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
        RuntimeAttr? runtime_attr_override

    }

    scatter (contig in contigs) {

        File bed_input = "~{bed_dir}/*~{contig}.bed.gz"
        File vcf_input = "~{vcf_dir}/*~{contig}.noheader.vcf.gz"
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
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_override
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
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}