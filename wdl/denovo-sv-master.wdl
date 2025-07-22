version 1.0

import "https://raw.githubusercontent.com/broadinstitute/gatk-sv/refs/tags/v1.0.4/wdl/AnnotateVcf.wdl" as Annotate

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
        Array[File] wes_denovo
        File wes_flipbook_responses
        File wes_file_paths_to_fix
        String release
        String denovo_docker
        RuntimeAttr? runtime_attr_wes_merge_to_annotate

        File contig_list
        File external_af_ref_bed
        File noncoding_bed
        File par_bed
        File ped_file
        File protein_coding_gtf
        Array[String] external_af_population
        String external_af_ref_prefix
        String sv_pipeline_docker
        String sv_base_mini_docker
        String gatk_docker

        RuntimeAttr? runtime_attr_wes_update_annotations
    }

    call denovo_wes_merge_to_annotate {
        input:
                wes_denovo = wes_denovo,
                wes_flipbook_responses = wes_flipbook_responses,
                wes_file_paths_to_fix = wes_file_paths_to_fix,
                release = release,
                denovo_docker = denovo_docker,
                runtime_attr_override = runtime_attr_wes_merge_to_annotate
    }

    call Annotate.AnnotateVcf as denovo_wes_annotate {
        input:
            contig_list = contig_list,
            gatk_docker = gatk_docker,
            prefix = release,
            sv_base_mini_docker = sv_base_mini_docker,
            sv_per_shard = 5000,
            sv_pipeline_docker = sv_pipeline_docker,
            vcf = denovo_wes_merge_to_annotate.vcf_to_annotate,
            external_af_population = external_af_population,
            external_af_ref_bed = external_af_ref_bed,
            external_af_ref_prefix = external_af_ref_prefix,
            noncoding_bed = noncoding_bed,
            par_bed = par_bed,
            ped_file = ped_file,
            protein_coding_gtf = protein_coding_gtf
    }

    call denovo_wes_update_annotations {
        input:
            wes_annotated_vcf = denovo_wes_annotate.annotated_vcf,
            bed_to_annotate = denovo_wes_merge_to_annotate.bed_to_annotate,
            release = release,
            denovo_docker = denovo_docker,
            runtime_attr_override = runtime_attr_wes_update_annotations
    }

    output{
        File wes_denovo_final = denovo_wes_update_annotations.denovo_wes_final
    }
}

task denovo_wes_merge_to_annotate {
    input{
        Array[File] wes_denovo
        File wes_flipbook_responses
        File wes_file_paths_to_fix
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

        Rscript /src/variant-interpretation/scripts/wes_denovo_merge.R \
            -c ~{sep="," wes_denovo} \
            -r ~{release} \
            -o .

        bedtools merge -d 15000000 \
            -i denovo_wes_for_merging-~{release}.bed \
            -c 4 -o collapse -delim ',' > denovo_wes_for_merging-~{release}.merged.bed

        Rscript /src/variant-interpretation/scripts/wes_denovo_post_merge_reformat.R \
            -d denovo_wes-~{release}.bed \
            -m denovo_wes_for_merging-~{release}.merged.bed \
            -f ~{wes_flipbook_responses} \
            -i ~{wes_file_paths_to_fix} \
            -r ~{release} \
            -o .

        ##Prepare to annoate
        Rscript /src/variant-interpretation/scripts/svs_add_bed2vcf.R \
            -i denovo_wes-~{release}.for_annotation.bed \
            -o denovo_wes-~{release}.for_annotation.vcf

        vcf-sort denovo_wes-~{release}.for_annotation.vcf | \
            bgzip -c > denovo_wes-~{release}.for_annotation.sorted.vcf.gz

        tabix -p vcf denovo_wes-~{release}.for_annotation.sorted.vcf.gz
    >>>

    output {
        File bed_to_annotate = "denovo_wes-~{release}.for_annotation.bed"
        File vcf_to_annotate = "denovo_wes-~{release}.for_annotation.sorted.vcf.gz"
        File vcf_idx_to_annotate = "denovo_wes-~{release}.for_annotation.sorted.vcf.gz"
    }
}

task denovo_wes_update_annotations {
    input{
        File wes_annotated_vcf
        File bed_to_annotate
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

        svtk vcf2bed -i ALL --include-filters ~{wes_annotated_vcf} - | \
            bgzip -c > denovo_wes-~{release}.annotated.bed.gz

        Rscript /src/variant-interpretation/scripts/wes_denovo_update_annotations.R \
            -s ~{bed_to_annotate} \
            -b denovo_wes-~{release}.annotated.bed.gz \
            -o denovo_wes-~{release}_final.bed
    >>>

    output {
        File denovo_wes_final = "denovo_wes-~{release}_final.bed"
    }
}

task denovo_wgs_merge_to_annotate {
    input{
        Array[File] wgs_denovo
        File wgs_flipbook_responses
        File mosaics_manual
        File mosaics_denovo_scores
        File tlocs_manual
        File ped_file
        File quality_control
        File genomic_disorders_manual
        File aneuploidies_manual
        File remove_svs
        File add_svs
        File baf_file
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

        Rscript /src/variant-interpretation/scripts/wgs_denovo_merge.R \
            -d ~{sep="," wgs_denovo} \
            -o denovo_wgs-~{release}.txt \
            -f ~{wgs_flipbook_responses} \
            -p ~{ped_file} \
            -q ~{quality_control}

        #run bedtools merge
        sort -k1,1 -k2,2n denovo_wgs-~{release}_ref_for_merge.txt > denovo_wgs-~{release}_ref_for_merge_sort.txt
        bedtools merge -d 500000 \
            -i denovo_wgs-~{release}_ref_for_merge_sort.txt -c 4 -o collapse \
            -delim ',' > denovo_wgs-~{release}_merged.txt

        #Remove overlapping mosaics
        head -n+1 ~{mosaics_manual} > mosaics.ovl.bed
        ##Keep all mosaics overlapping the callset to flag them
        bedtools intersect -wao -a ~{mosaics_manual} -b denovo_wgs-~{release}_merged.txt >> mosaics.ovl.bed

        #Remove overlapping translocations
        head -n+1 ~{tlocs_manual} > tlocs.ovl.bed
        tail -n+2 ~{tlocs_manual} > tlocs_noheader.bed
        ##Keep all tlocs overlapping the callset to flag them
        bedtools intersect -wao -a tlocs_noheader.bed -b denovo_wgs-~{release}_merged.txt >> tlocs.ovl.bed

        #Remove overlapping genomic disorders
        tail -n+2 ~{genomic_disorders_manual} > gds_noheader.bed
        ##Keep all gds overlapping the callset to flag them
        head -n+1 ~{genomic_disorders_manual} > gds.ovl.bed
        bedtools intersect -wao -a gds_noheader.bed -b denovo_wgs-~{release}_merged.txt >> gds.ovl.bed

        #Fix after running bedtools merge
        Rscript /src/variant-interpretation/scripts/wgs_denovo_post_merge_reformat.R \
            -d denovo_wgs-~{release}.txt \
            -e denovo_wgs-~{release}_merged.txt \
            -o denovo_wgs-~{release}.for_annotation.txt \
            -m mosaics.ovl.bed \
            -w ~{mosaics_denovo_scores} \
            -t tlocs.ovl.bed \
            -g gds.ovl.bed \
            -a ~{aneuploidies_manual} \
            -r ~{remove_svs} \
            -c ~{add_svs} \
            -p ~{ped_file} \
            -q ~{quality_control} \
            -b ~{baf_file}

        #Additionally reformat to re-annotate all calls
        Rscript /src/variant-interpretation/scripts/wgs_denovo_for_annotation.R \
            -i denovo_wgs-~{release}.for_annotation.txt \
            -o denovo_wgs-~{release}.for_annotation.bed

        ##Run BED to VCF script
        Rscript /src/variant-interpretation/scripts/svs_add_bed2vcf.R \
            -i denovo_wgs-~{release}.for_annotation.bed \
            -o denovo_wgs-~{release}.for_annotation.vcf

        vcf-sort denovo_wgs-~{release}.for_annotation.vcf > denovo_wgs-~{release}.for_annotation.sorted.vcf
        bgzip -c denovo_wgs-~{release}.for_annotation.sorted.vcf > denovo_wgs-~{release}.for_annotation.sorted.vcf.gz
        tabix -p vcf denovo_wgs-~{release}.for_annotation.sorted.vcf.gz
    >>>

    output {
        File bed_to_annotate = "denovo_wgs-~{release}.for_annotation.txt"
        File vcf_to_annotate = "denovo_wgs-~{release}.for_annotation.sorted.vcf.gz"
        File vcf_to_annotate_index = "denovo_wgs-~{release}.for_annotation.sorted.vcf.gz.tbi"
    }
}

task denovo_wgs_update_annotations {
    input{
        File wgs_annotated_vcf
        File bed_to_annotate
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
        ##Transform VCF annotated to BED file
        svtk vcf2bed --include-filters \
            -i ALL ~{wgs_annotated_vcf} - | \
            bgzip -c > denovo_wgs-~{release}.annotated.bed.gz

        ##Merge with master final bed file
        Rscript /src/variant-interpretation/scripts/denovosv_wgs_overwrite_predicted_csq.R \
            -s ~{bed_to_annotate} \
            -b denovo_wgs-~{release}.annotated.bed.gz \
            -o denovo_wgs-~{release}_final.bed
    >>>

    output {
        File denovo_wgs_final = "denovo_wgs-~{release}_final.bed"
    }
}