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
        String sv_base_mini_docker
        File hg38_fasta
        File hg38_fasta_fai
        File human_ancestor_fa
        File human_ancestor_fa_fai
        File top_level_fa
        File gerp_conservation_scores
        File hg38_vep_cache
        File loeuf_data
        String cohort_prefix
        RuntimeAttr? runtime_attr_normalize
        RuntimeAttr? runtime_attr_vep_annotate
        RuntimeAttr? runtime_attr_add_genotypes
    }

    call normalizeVCF {
        input:
            vcf_file=vcf_file,
            sv_base_mini_docker=sv_base_mini_docker,
            hg38_fasta=hg38_fasta,
            hg38_fasta_fai=hg38_fasta_fai,
            runtime_attr_override=runtime_attr_normalize
    }
        
    call vepAnnotate {
        input:
            vcf_file=normalizeVCF.vcf_no_genotype,
            vep_annotate_hail_python_script=vep_annotate_hail_python_script,
            top_level_fa=top_level_fa,
            human_ancestor_fa=human_ancestor_fa,
            human_ancestor_fa_fai=human_ancestor_fa_fai,
            gerp_conservation_scores=gerp_conservation_scores,
            hg38_vep_cache=hg38_vep_cache,
            loeuf_data=loeuf_data,
            vep_hail_docker=vep_hail_docker,
            runtime_attr_override=runtime_attr_vep_annotate
    }

    call addGenotypes { 
        input:
            vep_annotated_vcf=vepAnnotate.vep_vcf_file,
            normalized_vcf=normalizeVCF.vcf_normalized_file_with_genotype,
            normalized_vcf_idx=normalizeVCF.vcf_normalized_file_with_genotype_idx,
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_add_genotypes
    }

    output {   
        File vep_vcf_file = addGenotypes.merged_vcf_file
        File vep_vcf_idx = addGenotypes.merged_vcf_idx
    }
}   

task normalizeVCF {
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
    String vcf_normalized_file_name = "~{prefix}.normalized.vcf.gz"
    String vcf_normalized_nogeno_file_name = "~{prefix}.normalized.stripped.vcf.gz"

    command <<<
        set -euo pipefail

        bcftools norm -m - -f ~{hg38_fasta} -Oz -o ~{vcf_normalized_file_name} ~{vcf_file}
        bcftools index -t ~{vcf_normalized_file_name}

        bcftools view -G ~{vcf_normalized_file_name} -Oz -o ~{vcf_normalized_nogeno_file_name}
        bcftools index -t ~{vcf_normalized_nogeno_file_name}
        
    >>>
        
    output {
        File vcf_normalized_file_with_genotype = "~{prefix}.normalized.vcf.gz"
        File vcf_normalized_file_with_genotype_idx = "~{prefix}.normalized.vcf.gz.tbi"
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
        String vep_hail_docker
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(vcf_file, ".vcf.gz")
    String vep_annotated_vcf_name = "~{prefix}.vep.loftee.vcf.gz"

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
    String ancestor_dir = basename(human_ancestor_fa, "Homo_sapiens.GRCh38.dna.toplevel.fa.gz")

    command <<<
        set -euo pipefail

        dir_cache=$(dirname "~{hg38_vep_cache}")
        tar xzf ~{hg38_vep_cache} -C $dir_cache
        tabix -f -s 76 -b 77 -e 78 ~{loeuf_data}
        mv ~{gerp_conservation_scores} /opt/vep/Plugins/

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
        "--plugin", "LOEUF,file=~{loeuf_data},match_by=transcript",
        "--plugin", "LoF,loftee_path:/opt/vep/Plugins/,human_ancestor_fa:~{human_ancestor_fa},gerp_score:/opt/vep/Plugins/gerp_conservation_scores.homo_sapiens.GRCh38.bw",
        "-o", "STDOUT"],
        "vep_json_schema": "Struct{assembly_name:String,allele_string:String,ancestral:String,colocated_variants:Array[Struct{aa_allele:String,aa_maf:Float64,afr_allele:String,afr_maf:Float64,allele_string:String,amr_allele:String,amr_maf:Float64,clin_sig:Array[String],end:Int32,eas_allele:String,eas_maf:Float64,ea_allele:String,ea_maf:Float64,eur_allele:String,eur_maf:Float64,exac_adj_allele:String,exac_adj_maf:Float64,exac_allele:String,exac_afr_allele:String,exac_afr_maf:Float64,exac_amr_allele:String,exac_amr_maf:Float64,exac_eas_allele:String,exac_eas_maf:Float64,exac_fin_allele:String,exac_fin_maf:Float64,exac_maf:Float64,exac_nfe_allele:String,exac_nfe_maf:Float64,exac_oth_allele:String,exac_oth_maf:Float64,exac_sas_allele:String,exac_sas_maf:Float64,id:String,minor_allele:String,minor_allele_freq:Float64,phenotype_or_disease:Int32,pubmed:Array[Int32],sas_allele:String,sas_maf:Float64,somatic:Int32,start:Int32,strand:Int32}],context:String,end:Int32,id:String,input:String,intergenic_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],impact:String,minimised:Int32,variant_allele:String}],most_severe_consequence:String,motif_feature_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],high_inf_pos:String,impact:String,minimised:Int32,motif_feature_id:String,motif_name:String,motif_pos:Int32,motif_score_change:Float64,strand:Int32,variant_allele:String}],regulatory_feature_consequences:Array[Struct{allele_num:Int32,biotype:String,consequence_terms:Array[String],impact:String,minimised:Int32,regulatory_feature_id:String,variant_allele:String}],seq_region_name:String,start:Int32,strand:Int32,transcript_consequences:Array[Struct{allele_num:Int32,amino_acids:String,biotype:String,canonical:Int32,ccds:String,cdna_start:Int32,cdna_end:Int32,cds_end:Int32,cds_start:Int32,codons:String,consequence_terms:Array[String],distance:Int32,domains:Array[Struct{db:String,name:String}],exon:String,gene_id:String,gene_pheno:Int32,gene_symbol:String,gene_symbol_source:String,hgnc_id:String,hgvsc:String,hgvsp:String,hgvs_offset:Int32,impact:String,intron:String,lof:String,lof_flags:String,lof_filter:String,lof_info:String,minimised:Int32,polyphen_prediction:String,polyphen_score:Float64,protein_end:Int32,protein_start:Int32,protein_id:String,sift_prediction:String,sift_score:Float64,strand:Int32,swissprot:String,transcript_id:String,trembl:String,uniparc:String,variant_allele:String}],variant_class:String}"
        }' > vep_config.json

        python3.9 ~{vep_annotate_hail_python_script} ~{vcf_file} ~{vep_annotated_vcf_name} ~{cpu_cores} ~{memory}
        cp $(ls . | grep hail*.log) hail_log.txt
    >>>

    output {
        File vep_vcf_file = vep_annotated_vcf_name
        File hail_log = "hail_log.txt"
    }
}

task addGenotypes {
    input {
        File vep_annotated_vcf
        File normalized_vcf
        File normalized_vcf_idx
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
        Int? thread_num_override  
    }

    #  CleanVcf5.FindRedundantMultiallelics
    Float vep_annotate_sizes = size(vep_annotated_vcf, "GB") 
    Float norm_vcf_sizes = size(normalized_vcf, "GB")
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
        ~{normalized_vcf} \
        ~{vep_annotated_vcf}

        bcftools index -t ~{combined_vcf_name}

    >>>

    output {
        File merged_vcf_file = combined_vcf_name
        File merged_vcf_idx = combined_vcf_name + ".tbi"
    }
}