version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow DenovoSV_MASTER{
    input{
        Array[File] denovo_wes
        File flipbook_responses
        File file_paths_to_fix
        String release
        String denovo_docker
        RuntimeAttr? runtime_attr_wes_merge_to_annotate
    }
    call denovo_wes_merge_to_annotate {
        input:
                denovo_wes = denovo_wes,
                flipbook_responses = flipbook_responses,
                file_paths_to_fix = file_paths_to_fix,
                release = release,
                denovo_docker = denovo_docker,
                runtime_attr_override = runtime_attr_wes_merge_to_annotate
    }

    output{
        File wes_vcf_to_annotate = denovo_wes_merge_to_annotate.vcf_to_annotate
        File wes_vcf_idx_to_annotate = denovo_wes_merge_to_annotate.vcf_to_annotate
    }
}

task denovo_wes_merge_to_annotate {
    input{
        Array[File] denovo_wes
        File flipbook_responses
        File file_paths_to_fix
        String release
        String denovo_docker
        RuntimeAttr? runtime_attr_override
    }

    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    Float input_size = 12

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])

    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: denovo_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -eu

        Rscript wes_denovo_merge.R \
            -c ~{sep="," denovo_wes},
            -r ~{release} \
            -o .

        bedtools merge -d 15000000 \
            -i denovo_wes_for_merging-~{release}.bed \
            -c 4 -o collapse -delim ',' > denovo_wes_for_merging-~{release}.merged.bed

        Rscript wes_denovo_post_merge_reformat.R \
            -d denovo_wes-~{release}.bed \
            -m denovo_wes_for_merging-~{release}.merged.bed \
            -f ~{flipbook_responses} \
            -i ~{file_paths_to_fix} \
            -r ~{release} \
            -o .

        ##Prepare to annoate
        Rscript svs_add_bed2vcf.R \
            -i denovo_wes-~{release}.for_annotation.bed \
            -o denovo_wes-~{release}.for_annotation.vcf

        vcf-sort denovo_wes-~{release}.for_annotation.vcf | \
            bgzip -c > denovo_wes-~{release}.for_annotation.sorted.vcf.gz

        tabix -p vcf denovo_wes-~{release}.for_annotation.sorted.vcf.gz
    >>>

    output {
        File vcf_to_annotate = "denovo_wes-~{release}.for_annotation.sorted.vcf.gz"
        File vcf_idx_to_annotate = "denovo_wes-~{release}.for_annotation.sorted.vcf.gz"
    }
}


#$ svtk vcf2bed -i ALL --include-filters denovo_wes_merged-${release}.sorted.annotated.bed.gz
##update annotations
#$ Rscript denovosv_wgs_overwrite_predicted_csq.R -s denovo_wes_merged-${release}_final.bed
