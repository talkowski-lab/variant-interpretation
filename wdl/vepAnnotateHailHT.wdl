version 1.0
    
import "wes-denovo-helpers.wdl" as helpers

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
        String ht_uri
        String bucket_id

        File top_level_fa
        File ref_vep_cache

        File alpha_missense_file
        File eve_data

        String vep_hail_docker
        String genome_build='GRCh38'
        String vep_annotate_hail_ht_python_script

        # from extra
        String loeuf_v2_uri
        String loeuf_v4_uri
        File revel_file
        File clinvar_vcf_uri
        File omim_uri
        
        String gene_list='NA'
        String mpc_ht_uri
    }

    call helpers.getHailMTSize as getInputHTSize {
        input:
            mt_uri=ht_uri,
            hail_docker=vep_hail_docker
    }

    call vepAnnotate {
        input:
        ht_uri=ht_uri,
        bucket_id=bucket_id,
        input_size=getInputHTSize.mt_size,
        top_level_fa=top_level_fa,
        ref_vep_cache=ref_vep_cache,
        alpha_missense_file=alpha_missense_file,
        alpha_missense_file_idx=alpha_missense_file+'.tbi',
        eve_data=eve_data,
        eve_data_idx=eve_data+'.tbi',
        vep_hail_docker=vep_hail_docker,
        genome_build=genome_build,
        vep_annotate_hail_ht_python_script=vep_annotate_hail_ht_python_script,
        loeuf_v2_uri=loeuf_v2_uri,
        loeuf_v4_uri=loeuf_v4_uri,
        revel_file=revel_file,
        revel_file_idx=revel_file+'.tbi',
        clinvar_vcf_uri=clinvar_vcf_uri,
        omim_uri=omim_uri,
        gene_list=select_first([gene_list, 'NA']),
        mpc_ht_uri=mpc_ht_uri
    }

    output {
        String vep_annot_ht = vepAnnotate.vep_annot_ht
    }
}

task vepAnnotate {
    input {
        String ht_uri
        String bucket_id
        Float input_size

        File top_level_fa
        # File human_ancestor_fa
        # File human_ancestor_fa_fai
        # File gerp_conservation_scores
        File ref_vep_cache

        File alpha_missense_file
        File alpha_missense_file_idx
        File eve_data
        File eve_data_idx

        String vep_hail_docker
        String genome_build
        String vep_annotate_hail_ht_python_script

        # from extra
        String loeuf_v2_uri
        String loeuf_v4_uri
        File revel_file
        File revel_file_idx
        File clinvar_vcf_uri
        File omim_uri
        
        String gene_list
        String mpc_ht_uri

        RuntimeAttr? runtime_attr_override
    }

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

    command <<<
        set -euo pipefail

        dir_cache=$(dirname "~{ref_vep_cache}")
        tar xzf ~{ref_vep_cache} -C $dir_cache

        echo '{"command": [
        "/opt/vep/ensembl-vep/vep",
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
        "--assembly", "~{genome_build}",
        "--fasta", "~{top_level_fa}",
        "--plugin", "AlphaMissense,file=~{alpha_missense_file}",
        "--plugin", "EVE,file=~{eve_data}",
        "-o", "STDOUT"],
        "env": {
        "PERL5LIB": "/opt/vep/Plugins/"
        },
        "vep_json_schema": "Struct{assembly_name:String,allele_string:String,ancestral:String,colocated_variants:Array[Struct{aa_allele:String,aa_maf:Float64,afr_allele:String,afr_maf:Float64,allele_string:String,amr_allele:String,amr_maf:Float64,clin_sig:Array[String],end:Int32,eas_allele:String,eas_maf:Float64,ea_allele:String,ea_maf:Float64,eur_allele:String,eur_maf:Float64,exac_adj_allele:String,exac_adj_maf:Float64,exac_allele:String,exac_afr_allele:String,exac_afr_maf:Float64,exac_amr_allele:String,exac_amr_maf:Float64,exac_eas_allele:String,exac_eas_maf:Float64,exac_fin_allele:String,exac_fin_maf:Float64,exac_maf:Float64,exac_nfe_allele:String,exac_nfe_maf:Float64,exac_oth_allele:String,exac_oth_maf:Float64,exac_sas_allele:String,exac_sas_maf:Float64,id:String,minor_allele:String,minor_allele_freq:Float64,phenotype_or_disease:Int32,pubmed:Array[Int32],sas_allele:String,sas_maf:Float64,somatic:Int32,start:Int32,strand:Int32}],context:String,end:Int32,id:String,input:String,intergenic_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],impact:String,minimised:Int32,variant_allele:String}],most_severe_consequence:String,motif_feature_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],high_inf_pos:String,impact:String,minimised:Int32,motif_feature_id:String,motif_name:String,motif_pos:Int32,motif_score_change:Float64,strand:Int32,variant_allele:String}],regulatory_feature_consequences:Array[Struct{allele_num:Int32,biotype:String,consequence_terms:Array[String],impact:String,minimised:Int32,regulatory_feature_id:String,variant_allele:String}],seq_region_name:String,start:Int32,strand:Int32,transcript_consequences:Array[Struct{allele_num:Int32,amino_acids:String,biotype:String,canonical:Int32,ccds:String,cdna_start:Int32,cdna_end:Int32,cds_end:Int32,cds_start:Int32,codons:String,consequence_terms:Array[String],distance:Int32,domains:Array[Struct{db:String,name:String}],exon:String,gene_id:String,gene_pheno:Int32,gene_symbol:String,gene_symbol_source:String,hgnc_id:String,hgvsc:String,hgvsp:String,hgvs_offset:Int32,impact:String,intron:String,lof:String,lof_flags:String,lof_filter:String,lof_info:String,minimised:Int32,polyphen_prediction:String,polyphen_score:Float64,protein_end:Int32,protein_start:Int32,protein_id:String,sift_prediction:String,sift_score:Float64,strand:Int32,swissprot:String,transcript_id:String,trembl:String,uniparc:String,variant_allele:String}],variant_class:String}"
        }' > vep_config.json

        curl ~{vep_annotate_hail_ht_python_script} > vep_annotate.py
        proj_id=$(gcloud config get-value project)
        python3.9 vep_annotate.py -i ~{ht_uri} --bucket-id ~{bucket_id} --cores ~{cpu_cores} --mem ~{memory} \
        --build ~{genome_build} --project-id $proj_id --loeuf-v2 ~{loeuf_v2_uri} --loeuf-v4 ~{loeuf_v4_uri} \
        --mpc ~{mpc_ht_uri} --clinvar ~{clinvar_vcf_uri} --omim ~{omim_uri} \
        --revel ~{revel_file} --genes ~{gene_list} 
        cp $(ls . | grep hail*.log) hail_log.txt
    >>>

    output {
        File vep_annot_ht = read_lines('ht_uri.txt')[0]
        File hail_log = "hail_log.txt"
    }
}