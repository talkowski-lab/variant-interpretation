version 1.0
    
struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow vepAnnotateHail {

    input {
        File vcf_file
        File vep_annotate_hail_python_script
        String vep_hail_docker
        File hg38_fasta
        File hg38_fasta_fai
        File human_ancestor_fa
        File human_ancestor_fa_fai
        File top_level_fa
        File gerp_conservation_scores
        String cohort_prefix
        RuntimeAttr? runtime_attr_vep_annotate
    }

    call vepAnnotate {
        input:
        vcf_file=vcf_file,
        vep_annotate_hail_python_script=vep_annotate_hail_python_script,
        top_level_fa=top_level_fa,
        human_ancestor_fa=human_ancestor_fa,
        human_ancestor_fa_fai=human_ancestor_fa_fai,
        gerp_conservation_scores=gerp_conservation_scores,
        vep_hail_docker=vep_hail_docker,
        runtime_attr_override=runtime_attr_vep_annotate
    }

    output {   
        File vep_vcf_file = vepAnnotate.vep_vcf_file
    }
}   

task vepAnnotate {
    input {
        File vcf_file
        File vep_annotate_hail_python_script
        File top_level_fa
        File human_ancestor_fa
        File human_ancestor_fa_fai
        File gerp_conservation_scores
        String vep_hail_docker
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(vcf_file, ".vcf.gz")
    String vep_annotated_vcf_name = "~{prefix}.vep.loftee.vcf.gz"

    #  generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
    #  generally assume working memory is ~3 * inputs
    #  from task FinalCleanup in CleanVcfChromosome.wdl
    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    RuntimeAttr runtime_default = object {
        mem_gb: 8,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: vep_hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

    output {
        File vep_vcf_file = vep_annotated_vcf_name
    }

    String ancestor_dir = basename(human_ancestor_fa, "Homo_sapiens.GRCh38.dna.toplevel.fa.gz")

    command <<<
        set -euo pipefail

        echo '{"command": [
        "vep",
        "--format", "vcf",
        "__OUTPUT_FORMAT_FLAG__",
        "--force_overwrite",
        "-dir", "/opt/vep",
        "--everything",
        "--allele_number",
        "--no_stats",
        "--cache", "--offline",
        "--minimal",
        "--assembly", "GRCh38",
        "--fasta", "~{top_level_fa}",
        "--plugin", "LoF,loftee_path:/opt/vep/Plugins/,human_ancestor_fa:~{human_ancestor_fa},gerp_score:~{gerp_conservation_scores}",
        "--dir_plugins", "/opt/vep/Plugins/",
        "-o", "STDOUT"]
        }' > vep_config.json

        python3.9 ~{vep_annotate_hail_python_script} ~{vcf_file} ~{vep_annotated_vcf_name} 

    >>>
}
