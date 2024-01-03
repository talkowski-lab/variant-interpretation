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
        File hg38_fasta
        File hg38_fasta_fai
        File human_ancestor_fa
        File human_ancestor_fa_fai
        File top_level_fa
        File gerp_conservation_scores
        File hg38_vep_cache
        File loeuf_data
        File mpc_file
        String cohort_prefix
        String vep_hail_docker
        String sv_base_mini_docker
        Int? records_per_shard
        Int? thread_num_override
        RuntimeAttr? runtime_attr_split_vcf
        RuntimeAttr? runtime_attr_remove_genotypes
        RuntimeAttr? runtime_attr_vep_annotate
        RuntimeAttr? runtime_attr_add_genotypes
    }

    # call removeGenotypes {
    #     input:
    #         vcf_file=vcf_file,
    #         sv_base_mini_docker=sv_base_mini_docker,
    #         hg38_fasta=hg38_fasta,
    #         hg38_fasta_fai=hg38_fasta_fai,
    #         runtime_attr_override=runtime_attr_remove_genotypes
    # }
        
    # call vepAnnotate {
    #     input:
    #         vcf_file=removeGenotypes.vcf_no_genotype,
    #         vep_annotate_hail_python_script=vep_annotate_hail_python_script,
    #         top_level_fa=top_level_fa,
    #         human_ancestor_fa=human_ancestor_fa,
    #         human_ancestor_fa_fai=human_ancestor_fa_fai,
    #         gerp_conservation_scores=gerp_conservation_scores,
    #         hg38_vep_cache=hg38_vep_cache,
    #         loeuf_data=loeuf_data,
    #         mpc_file=mpc_file,
    #         vep_hail_docker=vep_hail_docker,
    #         runtime_attr_override=runtime_attr_vep_annotate
    # }

    # call addGenotypes { 
    #     input:
    #         vep_annotated_vcf=vepAnnotate.vep_vcf_file,
    #         vcf_file=vcf_file,
    #         vcf_file_idx=removeGenotypes.vcf_file_idx,
    #         sv_base_mini_docker=sv_base_mini_docker,
    #         runtime_attr_override=runtime_attr_add_genotypes
    # }

    # output {   
    #     File vep_vcf_file = addGenotypes.merged_vcf_file
    #     File vep_vcf_idx = addGenotypes.merged_vcf_idx
    # }
    
    if (defined(records_per_shard)) {
        call scatterVCF {
            input:
                vcf_file=vcf_file,
                records_per_shard=select_first([records_per_shard]),
                sv_base_mini_docker=sv_base_mini_docker,
                thread_num_override=thread_num_override,
                runtime_attr_override=runtime_attr_split_vcf
        }
    }
    
    if (!defined(records_per_shard)) {
        call splitByChromosome {
            input:
                vcf_file=vcf_file,
                sv_base_mini_docker=sv_base_mini_docker,
                thread_num_override=thread_num_override,
                runtime_attr_override=runtime_attr_split_vcf
        }
    }

    Array[File] vcf_shards = select_first([scatterVCF.shards, splitByChromosome.shards])

    scatter (vcf_shard in vcf_shards) {
        call vepAnnotate {
            input:
                vcf_file=vcf_shard,
                vep_annotate_hail_python_script=vep_annotate_hail_python_script,
                top_level_fa=top_level_fa,
                human_ancestor_fa=human_ancestor_fa,
                human_ancestor_fa_fai=human_ancestor_fa_fai,
                gerp_conservation_scores=gerp_conservation_scores,
                hg38_vep_cache=hg38_vep_cache,
                loeuf_data=loeuf_data,
                mpc_file=mpc_file,
                vep_hail_docker=vep_hail_docker,
                runtime_attr_override=runtime_attr_vep_annotate
        }
    }

    output {
        Array[File] vep_vcf_files = vepAnnotate.vep_vcf_file
        Array[File] vep_vcf_idx = vepAnnotate.vep_vcf_idx
        # File vep_vcf_file = vepAnnotate.vep_vcf_file
        # File vep_vcf_idx = vepAnnotate.vep_vcf_idx
    }
}   

task splitByChromosome { 
    input {
        File vcf_file
        String sv_base_mini_docker
        Int? thread_num_override
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    Int thread_num = select_first([thread_num_override,1])
    String prefix = basename(vcf_file, ".vcf.gz")

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
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

    command <<<
        chr_string="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
        echo $chr_string | tr ',' '\n' > chr_list.txt
        awk '$0="~{prefix}."$0".vcf.gz"' chr_list.txt > filenames.txt
        paste chr_list.txt filenames.txt > chr_filenames.txt
        bcftools +scatter ~{vcf_file} -Oz -o . --threads ~{thread_num} -S chr_filenames.txt 
    >>>

    output {
        Array[File] shards = glob("~{prefix}.chr*.vcf.gz")
        Array[String] shards_string = glob("~{prefix}.chr*.vcf.gz")
    }
}

task scatterVCF {
    input {
        File vcf_file
        Int records_per_shard
        String sv_base_mini_docker
        Int? thread_num_override
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    Int thread_num = select_first([thread_num_override,1])
    String prefix = basename(vcf_file, ".vcf.gz")

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
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
    
    command <<<
        set -euo pipefail
        #  in case the file is empty create an empty shard
        bcftools view -h ~{vcf_file} | bgzip -c > "~{prefix}.0.vcf.gz"
        bcftools +scatter ~{vcf_file} -o . -Oz -p "~{prefix}". --threads ~{thread_num} -n ~{records_per_shard}

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

task removeGenotypes {
    input {
        File vcf_file
        String sv_base_mini_docker
        File hg38_fasta
        File hg38_fasta_fai
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    RuntimeAttr runtime_default = object {
        mem_gb: 4,
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
    String vcf_normalized_nogeno_file_name = "~{prefix}.normalized.stripped.vcf.gz"

    command <<<
        set -euo pipefail

        bcftools index -t ~{vcf_file} -o ~{basename(vcf_file)+'.tbi'}
        bcftools view -G ~{vcf_file} -Oz -o ~{vcf_normalized_nogeno_file_name}
        bcftools index -t ~{vcf_normalized_nogeno_file_name}
        
    >>>
        
    output {
        File vcf_file_idx = basename(vcf_file) + '.tbi'
        File vcf_no_genotype = "~{prefix}.normalized.stripped.vcf.gz"
        File vcf_no_genotype_idx = "~{prefix}.normalized.stripped.vcf.gz.tbi"
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
        File hg38_vep_cache
        File loeuf_data
        File mpc_file
        String vep_hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB") + size(hg38_vep_cache, "GB") + size(gerp_conservation_scores, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 10.0
    RuntimeAttr runtime_default = object {
        mem_gb: 8,
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
        docker: vep_hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])
    
    String prefix = basename(vcf_file, ".vcf.gz")
    String vep_annotated_vcf_name = "~{prefix}.vep.vcf.bgz"

    command <<<
        set -euo pipefail

        dir_cache=$(dirname "~{hg38_vep_cache}")
        tar xzf ~{hg38_vep_cache} -C $dir_cache
        tabix -f -s 76 -b 77 -e 78 ~{loeuf_data}

        echo '{"command": [
        "vep",
        "--format", "vcf",
        "__OUTPUT_FORMAT_FLAG__",
        "--force_overwrite",
        "--dir_plugins", "/opt/vep/Plugins/",
        "-dir_cache", "'$dir_cache'",
        "--everything",
        "--allele_number",
        "--no_stats",
        "--cache", 
        "--offline",
        "--minimal",
        "--assembly", "GRCh38",
        "--fasta", "~{top_level_fa}",
        "--plugin", "MPC,~{mpc_file}",
        "--plugin", "LOEUF,file=~{loeuf_data},match_by=transcript",
        "--plugin", "LoF,loftee_path:/opt/vep/Plugins/,human_ancestor_fa:~{human_ancestor_fa},gerp_bigwig:~{gerp_conservation_scores}",
        "-o", "STDOUT"],
        "vep_json_schema": "Struct{assembly_name:String,allele_string:String,ancestral:String,colocated_variants:Array[Struct{aa_allele:String,aa_maf:Float64,afr_allele:String,afr_maf:Float64,allele_string:String,amr_allele:String,amr_maf:Float64,clin_sig:Array[String],end:Int32,eas_allele:String,eas_maf:Float64,ea_allele:String,ea_maf:Float64,eur_allele:String,eur_maf:Float64,exac_adj_allele:String,exac_adj_maf:Float64,exac_allele:String,exac_afr_allele:String,exac_afr_maf:Float64,exac_amr_allele:String,exac_amr_maf:Float64,exac_eas_allele:String,exac_eas_maf:Float64,exac_fin_allele:String,exac_fin_maf:Float64,exac_maf:Float64,exac_nfe_allele:String,exac_nfe_maf:Float64,exac_oth_allele:String,exac_oth_maf:Float64,exac_sas_allele:String,exac_sas_maf:Float64,id:String,minor_allele:String,minor_allele_freq:Float64,phenotype_or_disease:Int32,pubmed:Array[Int32],sas_allele:String,sas_maf:Float64,somatic:Int32,start:Int32,strand:Int32}],context:String,end:Int32,id:String,input:String,intergenic_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],impact:String,minimised:Int32,variant_allele:String}],most_severe_consequence:String,motif_feature_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],high_inf_pos:String,impact:String,minimised:Int32,motif_feature_id:String,motif_name:String,motif_pos:Int32,motif_score_change:Float64,strand:Int32,variant_allele:String}],regulatory_feature_consequences:Array[Struct{allele_num:Int32,biotype:String,consequence_terms:Array[String],impact:String,minimised:Int32,regulatory_feature_id:String,variant_allele:String}],seq_region_name:String,start:Int32,strand:Int32,transcript_consequences:Array[Struct{allele_num:Int32,amino_acids:String,biotype:String,canonical:Int32,ccds:String,cdna_start:Int32,cdna_end:Int32,cds_end:Int32,cds_start:Int32,codons:String,consequence_terms:Array[String],distance:Int32,domains:Array[Struct{db:String,name:String}],exon:String,gene_id:String,gene_pheno:Int32,gene_symbol:String,gene_symbol_source:String,hgnc_id:String,hgvsc:String,hgvsp:String,hgvs_offset:Int32,impact:String,intron:String,lof:String,lof_flags:String,lof_filter:String,lof_info:String,minimised:Int32,polyphen_prediction:String,polyphen_score:Float64,protein_end:Int32,protein_start:Int32,protein_id:String,sift_prediction:String,sift_score:Float64,strand:Int32,swissprot:String,transcript_id:String,trembl:String,uniparc:String,variant_allele:String}],variant_class:String}"
        }' > vep_config.json

        python3.9 ~{vep_annotate_hail_python_script} ~{vcf_file} ~{vep_annotated_vcf_name} ~{cpu_cores} ~{memory}
        cp $(ls . | grep hail*.log) hail_log.txt
        /opt/vep/bcftools/bcftools index -t ~{vep_annotated_vcf_name}
    >>>

    output {
        File vep_vcf_file = vep_annotated_vcf_name
        File vep_vcf_idx = vep_annotated_vcf_name + '.tbi'
        File hail_log = "hail_log.txt"
    }
}

task addGenotypes {
    input {
        File vep_annotated_vcf
        File vcf_file
        File vcf_file_idx
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
        Int? thread_num_override  
    }

    #  CleanVcf5.FindRedundantMultiallelics
    Float vep_annotate_sizes = size(vep_annotated_vcf, "GB") 
    Float norm_vcf_sizes = size(vcf_file, "GB")
    Float base_disk_gb = 10.0

    RuntimeAttr runtime_default = object {
                                      mem_gb: 4,  
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
        ~{vcf_file} \
        ~{vep_annotated_vcf}

        bcftools index -t ~{combined_vcf_name}

    >>>

    output {
        File merged_vcf_file = combined_vcf_name
        File merged_vcf_idx = combined_vcf_name + ".tbi"
    }
}