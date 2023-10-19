version 1.0
    
import "vepAnnotateSingle.wdl" as vepAnnotateSingle

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow vepAnnotate {
    input {
        # file can be a list of vcf files or just one vcf file
        File file
        String vep_docker
        String sv_base_mini_docker
        File hg38_fasta
        File hg38_fasta_fai
        File human_ancestor_fa
        File human_ancestor_fa_fai
        File top_level_fa
        Boolean merge_annotated_vcfs=true
        Int? records_per_shard
        RuntimeAttr? runtime_attr_normalize
        RuntimeAttr? runtime_attr_split_vcf
        RuntimeAttr? runtime_attr_vep_annotate
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_add_genotypes
    }

    String filename = basename(file)

    # if file is vcf.gz (just one file)
    if (sub(filename, ".vcf.gz", "") != filename) {
        Array[String] vcf_files = [file]    
    }
    if (sub(filename, ".vcf.gz", "") == filename) {
        Array[String] vcf_files = read_lines(file)
    }

    scatter (vcf_file in vcf_files) {
        File vcf_file = vcf_file
        call vepAnnotateSingle {
            input:
                vcf_file=vcf_file,
                vep_docker=vep_docker,
                sv_base_mini_docker=sv_base_mini_docker,
                hg38_fasta=hg38_fasta,
                hg38_fasta_fai=hg38_fasta_fai,
                human_ancestor_fa=human_ancestor_fa,
                human_ancestor_fa_fai=human_ancestor_fa_fai,
                top_level_fa=top_level_fa, 
                cohort_prefix=basename(vcf_file, ".vcf.gz"),
                merge_annotated_vcfs=merge_annotated_vcfs,
                records_per_shard=select_first([records_per_shard]),
                runtime_attr_normalize=runtime_attr_normalize,
                runtime_attr_split_vcf=runtime_attr_split_vcf,
                runtime_attr_vep_annotate=runtime_attr_vep_annotate,
                runtime_attr_merge_vcfs=runtime_attr_merge_vcfs,
                runtime_attr_add_genotypes=runtime_attr_add_genotypes
        }
    }

    output {
        Array[Array[File]] vep_annotated_final_vcf = vepAnnotateSingle.vep_annotated_final_vcf
        Array[Array[File]] vep_annotated_final_vcf_idx = vepAnnotateSingle.vep_annotated_final_vcf_idx
    }
}