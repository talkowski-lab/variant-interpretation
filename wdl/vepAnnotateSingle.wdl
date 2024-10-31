version 1.0
    
struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow vepAnnotateSingle {

    input {
        File vcf_file
        String vep_hail_docker
        String sv_base_mini_docker
        File hg38_fasta
        File hg38_fasta_fai
        File human_ancestor_fa
        File human_ancestor_fa_fai
        File top_level_fa
        File gerp_conservation_scores
        String cohort_prefix
        Boolean merge_annotated_vcfs
        Int? records_per_shard
        RuntimeAttr? runtime_attr_normalize
        RuntimeAttr? runtime_attr_split_vcf
        RuntimeAttr? runtime_attr_vep_annotate
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_add_genotypes
    }

    # split by shard, vep annotate with Loftee 
    if (defined(records_per_shard)) {
        
        call scatterVCF {
            input:
                vcf_file=vcf_file,
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override = runtime_attr_split_vcf,
                records_per_shard=select_first([records_per_shard]),
                sv_base_mini_docker=sv_base_mini_docker
        }

        scatter (shard in scatterVCF.shards) {
            # normalize vcf file
            call normalizeVCF as normalizeVCF_sharded {
                input:
                    vcf_file=shard,
                    sv_base_mini_docker=sv_base_mini_docker,
                    hg38_fasta=hg38_fasta,
                    hg38_fasta_fai=hg38_fasta_fai,
                    runtime_attr_override=runtime_attr_normalize
            }

            call vepAnnotate as vepAnnotate_sharded {
                input:
                    vcf_file=normalizeVCF_sharded.vcf_no_genotype,
                    top_level_fa=top_level_fa,
                    human_ancestor_fa=human_ancestor_fa,
                    human_ancestor_fa_fai=human_ancestor_fa_fai,
                    gerp_conservation_scores=gerp_conservation_scores,
                    vep_hail_docker=vep_hail_docker,
                    runtime_attr_override=runtime_attr_vep_annotate
            }

            call addGenotypes as addGenotypes_sharded { 
                input:
                vep_annotated_vcf=select_first([vepAnnotate_sharded.vep_vcf_file]),
                normalized_vcf=normalizeVCF_sharded.vcf_normalized_file_with_genotype,
                normalized_vcf_idx=normalizeVCF_sharded.vcf_normalized_file_with_genotype_idx,
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_add_genotypes
            }
        }
    }

    if (!defined(records_per_shard)) {
        # call normalizeVCF {
        #     input:
        #         vcf_file=vcf_file,
        #         sv_base_mini_docker=sv_base_mini_docker,
        #         hg38_fasta=hg38_fasta,
        #         hg38_fasta_fai=hg38_fasta_fai,
        #         runtime_attr_override=runtime_attr_normalize
        # }
        call vepAnnotate {
            input:
                # vcf_file=normalizeVCF.vcf_no_genotype,
                vcf_file=vcf_file,
                top_level_fa=top_level_fa,
                human_ancestor_fa=human_ancestor_fa,
                human_ancestor_fa_fai=human_ancestor_fa_fai,
                gerp_conservation_scores=gerp_conservation_scores,
                vep_hail_docker=vep_hail_docker,
                runtime_attr_override=runtime_attr_vep_annotate
        }
        # call addGenotypes { 
        #     input:
        #     vep_annotated_vcf=select_first([vepAnnotate.vep_vcf_file]),
        #     normalized_vcf=normalizeVCF.vcf_normalized_file_with_genotype,
        #     normalized_vcf_idx=normalizeVCF.vcf_normalized_file_with_genotype_idx,
        #     sv_base_mini_docker=sv_base_mini_docker,
        #     runtime_attr_override=runtime_attr_add_genotypes
        # }
    }

    # Array[File] genotyped_vcf_file = select_first([addGenotypes_sharded.merged_vcf_file, addGenotypes.merged_vcf_file])
    # Array[File] genotyped_vcf_idx = select_first([addGenotypes_sharded.merged_vcf_idx, addGenotypes.merged_vcf_idx])


    Array[File] genotyped_vcf_file = select_first([addGenotypes_sharded.merged_vcf_file, vepAnnotate.vep_vcf_file])
    Array[File] genotyped_vcf_idx = select_first([addGenotypes_sharded.merged_vcf_idx, vepAnnotate.vep_vcf_idx])

    if (merge_annotated_vcfs) {
        call mergeVCFs {
            input:
                vcf_contigs=genotyped_vcf_file,
                sv_base_mini_docker=sv_base_mini_docker,
                cohort_prefix=cohort_prefix,
                runtime_attr_override=runtime_attr_merge_vcfs
        }       
        File merged_vcf_file = select_first([mergeVCFs.merged_vcf_file])
        File merged_vcf_idx = select_first([mergeVCFs.merged_vcf_idx])
    }

    output {   
        Array[File] vep_annotated_final_vcf = select_first([merged_vcf_file, genotyped_vcf_file])
        Array[File] vep_annotated_final_vcf_idx = select_first([merged_vcf_idx, genotyped_vcf_idx])
    }
}   

task addGenotypes {
    input {
        File vep_annotated_vcf
        File normalized_vcf
        File normalized_vcf_idx
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
        Int? thread_num_override  # TODO: 
    }

    #  CleanVcf5.FindRedundantMultiallelics
    Float vep_annotate_sizes = size(vep_annotated_vcf, "GB") 
    Float norm_vcf_sizes = size(normalized_vcf, "GB")
    Float base_disk_gb = 10.0

    RuntimeAttr runtime_default = object {
                                      mem_gb: 16,  # TODO: lower
                                      disk_gb: ceil(base_disk_gb + (vep_annotate_sizes + norm_vcf_sizes) * 5.0),
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(vep_annotated_vcf, ".vcf.gz")
    String combined_vcf_name="~{prefix}.vep.geno.vcf.gz"
    Int thread_num = select_first([thread_num_override,1])

    command <<<
        set -euo pipefail
        
        bcftools index -t ~{vep_annotated_vcf}

        bcftools merge \
        --no-version \
        --threads ~{thread_num} \
        -Oz \
        --output ~{combined_vcf_name} \
        ~{normalized_vcf} \
        ~{vep_annotated_vcf}

        bcftools index -t ~{combined_vcf_name}

    >>>

    output {
        File merged_vcf_file = combined_vcf_name
        File merged_vcf_idx = combined_vcf_name + ".tbi"
    }

}


task mergeVCFs {
    input {
        Array[File] vcf_contigs
        String sv_base_mini_docker
        String cohort_prefix
        RuntimeAttr? runtime_attr_override
    }

    #  generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
    #  generally assume working memory is ~3 * inputs
    #  CleanVcf5.FindRedundantMultiallelics
    Float input_size = size(vcf_contigs, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 2.0
    Float input_mem_scale = 3.0
    Float input_disk_scale = 5.0  # TODO: change to 10?
    
    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb + input_size * input_mem_scale,  # TODO: fix, 4 GB
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"  # TODO: SSD
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String merged_vcf_name="~{cohort_prefix}.vep.merged.vcf.gz"

    command <<<
        set -euo pipefail
        VCFS="~{write_lines(vcf_contigs)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        bcftools concat --no-version --naive -Oz --file-list vcfs_sorted.list --output ~{merged_vcf_name}
        bcftools index -t ~{merged_vcf_name}
    >>>

    output {
        File merged_vcf_file=merged_vcf_name
        File merged_vcf_idx=merged_vcf_name + ".tbi"
    }
}

task normalizeVCF{
    input {
        File vcf_file
        String sv_base_mini_docker
        File hg38_fasta
        File hg38_fasta_fai
        RuntimeAttr? runtime_attr_override
    }

    #  generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
    #  generally assume working memory is ~3 * inputs
    #  from task FinalCleanup in CleanVcfChromosome.wdl
    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 2.0
    Float input_mem_scale = 3.0
    Float input_disk_scale = 5.0
    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb + input_size * input_mem_scale,
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(vcf_file, ".vcf.gz")
    String vcf_normalized_file_name = "~{prefix}.normalized.vcf.gz"
    String vcf_normalized_nogeno_file_name = "~{prefix}.normalized.stripped.vcf.gz"
    String vcf_normalized_nogeno_idx_name = "~{prefix}.normalized.stripped.vcf.gz.tbi"


    output{
        File vcf_normalized_file_with_genotype = "~{prefix}.normalized.vcf.gz"
        File vcf_normalized_file_with_genotype_idx = "~{prefix}.normalized.vcf.gz.tbi"
        File vcf_no_genotype = "~{prefix}.normalized.stripped.vcf.gz"
        File vcf_no_genotype_idx = "~{prefix}.normalized.stripped.vcf.gz.tbi"
    }

    command <<<
        set -euo pipefail

        bcftools norm -m - -f ~{hg38_fasta} -Oz -o ~{vcf_normalized_file_name} ~{vcf_file}

        bcftools index -t ~{vcf_normalized_file_name}

        bcftools view -G ~{vcf_normalized_file_name} -Oz -o ~{vcf_normalized_nogeno_file_name}

        bcftools index -t ~{vcf_normalized_nogeno_file_name}
        
    >>>

}

task scatterVCF {
    input {
        File vcf_file
        Int records_per_shard
        String? contig
        String sv_base_mini_docker
        Int? thread_num_override
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 2.0
    Float input_mem_scale = 3.0
    Float input_disk_scale = 5.0
    Int thread_num = select_first([thread_num_override,1])
    String prefix = basename(vcf_file, ".vcf.gz")

    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb + input_size * input_mem_scale,
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

    command <<<
        set -euo pipefail
        #  in case the file is empty create an empty shard
        bcftools view -h ~{vcf_file} | bgzip -c > "~{prefix}.0.vcf.gz"
        bcftools +scatter ~{vcf_file} -o . -O z -p "~{prefix}". --threads ~{thread_num} -n ~{records_per_shard} ~{"-r " + contig}

        ls "~{prefix}".*.vcf.gz | sort -k1,1V > vcfs.list
        i=0
        while read VCF; do
          shard_no=`printf %06d $i`
          mv "$VCF" "~{prefix}.shard_${shard_no}.vcf.gz"
          i=$((i+1))
        done < vcfs.list
    >>>
    output {
        Array[File] shards = glob("~{prefix}.shard_*.vcf.gz")
        Array[String] shards_string = glob("~{prefix}.shard_*.vcf.gz")
    }
}

task vepAnnotate {
    input {
        File vcf_file
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
    Float input_size = size(vcf_file, "GB") + size(hg38_vep_cache, "GB") + size(gerp_conservation_scores, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 2.0
    Float input_mem_scale = 3.0
    Float input_disk_scale = 5.0
    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb + input_size * input_mem_scale,
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
        File vep_vcf_idx = vep_annotated_vcf_name + '.tbi'
    }

    String ancestor_dir = basename(human_ancestor_fa, "Homo_sapiens.GRCh38.dna.toplevel.fa.gz")

    command <<<
        set -euo pipefail
        vep --vcf \
        --verbose \
        --force_overwrite \
        -dir /opt/vep/.vep \
        --format vcf \
        --everything \
        --allele_number \
        --no_stats \
        --dir_cache . \
        --cache \
        --offline \
        --minimal \
        --assembly GRCh38 \
        --fasta ~{top_level_fa} \
        --input_file ~{vcf_file} \
        --output_file ~{vep_annotated_vcf_name} \
        --compress_output bgzip \
        --plugin LoF,loftee_path:/opt/vep/Plugins/,human_ancestor_fa:~{human_ancestor_fa},gerp_bigwig:~{gerp_conservation_scores} \
        --dir_plugins /opt/vep/Plugins/

        /opt/vep/bcftools/bcftools index -t ~{vep_annotated_vcf_name}
    >>>
}
