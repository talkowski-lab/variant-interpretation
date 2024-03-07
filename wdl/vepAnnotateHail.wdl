version 1.0
    
import "scatterVCF.wdl" as scatterVCF
import "mergeSplitVCF.wdl" as mergeSplitVCF
import "mergeVCFs.wdl" as mergeVCFs

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
        File? vcf_file
        String? mt_uri
        String vep_annotate_hail_python_script
        String split_vcf_hail_script
        File hg38_fasta
        File hg38_fasta_fai
        File human_ancestor_fa
        File human_ancestor_fa_fai
        File top_level_fa
        File gerp_conservation_scores
        File hg38_vep_cache
        File loeuf_data
        String cohort_prefix
        String hail_docker
        String vep_hail_docker
        String sv_base_mini_docker
        Boolean split_by_chromosome
        Boolean split_into_shards 
        Boolean merge_split_vcf
        Int shards_per_chunk=10  # combine pre-sharded VCFs
        Array[File]? vcf_shards  # if scatterVCF.wdl already run before VEP
        Array[String]? row_fields_to_keep=[]
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_vep_annotate
        ## parameters for scatterVCF
        # Boolean has_index=false
        # Boolean localize_vcf=false
        # Int mt_size=0
        # Int n_shards=0
        # Int records_per_shard=0
    }

    String file = select_first([vcf_file, mt_uri])
    # input is not a single VCF file, so merge shards in chunks, then run VEP on merged chunks
    if (merge_split_vcf) { 
        # combine pre-sharded VCFs into chunks
        call mergeSplitVCF.splitFile as splitFile {
            input:
                file=file,
                shards_per_chunk=shards_per_chunk,
                cohort_prefix=cohort_prefix,
                vep_hail_docker=vep_hail_docker
        }
        scatter (chunk_file in splitFile.chunks) {
            call mergeVCFs.mergeVCFs as mergeVCFs {
                input:
                    vcf_files=read_lines(chunk_file),
                    sv_base_mini_docker=sv_base_mini_docker,
                    cohort_prefix=basename(chunk_file),
                    runtime_attr_override=runtime_attr_merge_vcfs
            }
            call vepAnnotate as vepAnnotateMergedShards {
                input:
                    vcf_file=mergeVCFs.merged_vcf_file,
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
        }
    }

    if (!merge_split_vcf) {
        if (!defined(vcf_shards)) {
            call scatterVCF.scatterVCF_workflow as scatterVCF {
                input:
                    file=file,
                    split_vcf_hail_script=split_vcf_hail_script,
                    cohort_prefix=cohort_prefix,
                    hail_docker=hail_docker,
                    sv_base_mini_docker=sv_base_mini_docker,
                    split_by_chromosome=split_by_chromosome,
                    split_into_shards=split_into_shards,
                    row_fields_to_keep=row_fields_to_keep
            }
        }
        Array[File] vcf_shards_ = select_first([scatterVCF.vcf_shards, vcf_shards])
    
        scatter (vcf_shard in vcf_shards_) {
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
                    vep_hail_docker=vep_hail_docker,
                    runtime_attr_override=runtime_attr_vep_annotate
            }
        }
    }

    Array[File] vep_vcf_files_ = select_first([vepAnnotate.vep_vcf_file, vepAnnotateMergedShards.vep_vcf_file])
    Array[File] vep_vcf_idx_ = select_first([vepAnnotate.vep_vcf_idx, vepAnnotateMergedShards.vep_vcf_idx])

    output {
        Array[File] vep_vcf_files = vep_vcf_files_
        Array[File] vep_vcf_idx = vep_vcf_idx_
    }
}   

task vepAnnotate {
    input {
        File vcf_file
        String vep_annotate_hail_python_script
        File top_level_fa
        File human_ancestor_fa
        File human_ancestor_fa_fai
        File gerp_conservation_scores
        File hg38_vep_cache
        File loeuf_data
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

    String filename = basename(vcf_file)
    String prefix = if (sub(filename, "\\.gz", "")!=filename) then basename(vcf_file, ".vcf.gz") else basename(vcf_file, ".vcf.bgz")
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
        "--plugin", "LOEUF,file=~{loeuf_data},match_by=transcript",
        "--plugin", "LoF,loftee_path:/opt/vep/Plugins/,human_ancestor_fa:~{human_ancestor_fa},gerp_bigwig:~{gerp_conservation_scores}",
        "-o", "STDOUT"],
        "env": {
        "PERL5LIB": "/opt/vep/Plugins/"
        },
        "vep_json_schema": "Struct{assembly_name:String,allele_string:String,ancestral:String,colocated_variants:Array[Struct{aa_allele:String,aa_maf:Float64,afr_allele:String,afr_maf:Float64,allele_string:String,amr_allele:String,amr_maf:Float64,clin_sig:Array[String],end:Int32,eas_allele:String,eas_maf:Float64,ea_allele:String,ea_maf:Float64,eur_allele:String,eur_maf:Float64,exac_adj_allele:String,exac_adj_maf:Float64,exac_allele:String,exac_afr_allele:String,exac_afr_maf:Float64,exac_amr_allele:String,exac_amr_maf:Float64,exac_eas_allele:String,exac_eas_maf:Float64,exac_fin_allele:String,exac_fin_maf:Float64,exac_maf:Float64,exac_nfe_allele:String,exac_nfe_maf:Float64,exac_oth_allele:String,exac_oth_maf:Float64,exac_sas_allele:String,exac_sas_maf:Float64,id:String,minor_allele:String,minor_allele_freq:Float64,phenotype_or_disease:Int32,pubmed:Array[Int32],sas_allele:String,sas_maf:Float64,somatic:Int32,start:Int32,strand:Int32}],context:String,end:Int32,id:String,input:String,intergenic_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],impact:String,minimised:Int32,variant_allele:String}],most_severe_consequence:String,motif_feature_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],high_inf_pos:String,impact:String,minimised:Int32,motif_feature_id:String,motif_name:String,motif_pos:Int32,motif_score_change:Float64,strand:Int32,variant_allele:String}],regulatory_feature_consequences:Array[Struct{allele_num:Int32,biotype:String,consequence_terms:Array[String],impact:String,minimised:Int32,regulatory_feature_id:String,variant_allele:String}],seq_region_name:String,start:Int32,strand:Int32,transcript_consequences:Array[Struct{allele_num:Int32,amino_acids:String,biotype:String,canonical:Int32,ccds:String,cdna_start:Int32,cdna_end:Int32,cds_end:Int32,cds_start:Int32,codons:String,consequence_terms:Array[String],distance:Int32,domains:Array[Struct{db:String,name:String}],exon:String,gene_id:String,gene_pheno:Int32,gene_symbol:String,gene_symbol_source:String,hgnc_id:String,hgvsc:String,hgvsp:String,hgvs_offset:Int32,impact:String,intron:String,lof:String,lof_flags:String,lof_filter:String,lof_info:String,minimised:Int32,polyphen_prediction:String,polyphen_score:Float64,protein_end:Int32,protein_start:Int32,protein_id:String,sift_prediction:String,sift_score:Float64,strand:Int32,swissprot:String,transcript_id:String,trembl:String,uniparc:String,variant_allele:String}],variant_class:String}"
        }' > vep_config.json

        curl ~{vep_annotate_hail_python_script} > vep_annotate.py
        python3.9 vep_annotate.py ~{vcf_file} ~{vep_annotated_vcf_name} ~{cpu_cores} ~{memory}
        cp $(ls . | grep hail*.log) hail_log.txt
        /opt/vep/bcftools/bcftools index -t ~{vep_annotated_vcf_name}
    >>>

    output {
        File vep_vcf_file = vep_annotated_vcf_name
        File vep_vcf_idx = vep_annotated_vcf_name + '.tbi'
        File hail_log = "hail_log.txt"
    }
}