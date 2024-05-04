##Create workflow from firecloud to add disable dict validation option

workflow CC_Crams{
    Array[File] crams
    Array[File] crais
    Array[String] sample_IDs
    File hg38_reference
    File hg38_reference_fai
    File hg38_reference_dict

    File editRGSMjar
    File intervals_exons_all
    File? gatk4_jar_override
    String gatk_docker

    scatter (scatter_index in range(length(crams))){
        call CollectCountsCram as coverage_all{
            input:
                intervals_exons = intervals_exons_all,
                cram = crams[scatter_index],
                crai = crais[scatter_index],
                sample_ID = sample_IDs[scatter_index],
                hg38_reference = hg38_reference,
                hg38_reference_fai = hg38_reference_fai,
                hg38_reference_dict = hg38_reference_dict,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                editRGSMjar = editRGSMjar
        }
    }
    output {
        Array [String] sample_all = coverage_all.entity_id
        Array [File] counts_exons_all = coverage_all.counts_exons
        Array [String] source_crams_all = coverage_all.source_cram
    }
}

task CollectCountsCram {
    File intervals_exons
    File cram
    File crai
    File hg38_reference
    File hg38_reference_fai
    File hg38_reference_dict
    File? gatk4_jar_override
    String sample_ID
    File editRGSMjar

    # Runtime parameters
    String gatk_docker
    Int? disk_space_gb
    Boolean use_ssd = false
    Int cpu=1
    Int? preemptible_attempts

    Int machine_mem_mb = 600
    Int command_mem_mb = machine_mem_mb - 100

    String base_filename = basename(cram, ".cram")
    String counts_exons_filename = "${sample_ID}.exons.counts.tsv"

    command <<<
        set -euo pipefail
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx${command_mem_mb}m" CollectReadCounts \
            -I ${cram} \
            --read-index ${crai} \
            -L ${intervals_exons} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --reference ${hg38_reference} \
            --format TSV \
            --disable-sequence-dictionary-validation true \
            -O temp_file.tsv

        java "-Xmx${command_mem_mb}m" -jar ${editRGSMjar} temp_file.tsv ${counts_exons_filename} ${sample_ID}

        rm -rf temp_file.tsv

    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(cram, "GB")) + 50]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
        maxRetries: 3
    }

    output {
        String entity_id = sample_ID
        File counts_exons = counts_exons_filename
        String source_cram = cram

    }
}