version 1.0
    
import "scatterVCF.wdl" as scatterVCF
import "mergeSplitVCF.wdl" as mergeSplitVCF
import "mergeVCFs.wdl" as mergeVCFs
import "vepAnnotateHailExtra.wdl" as annotateExtra

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

        File ref_fasta
        File ref_fasta_fai
        File human_ancestor_fa
        File human_ancestor_fa_fai
        File top_level_fa
        File gerp_conservation_scores
        File ref_vep_cache

        File alpha_missense_file
        File eve_data

        # File revel_file
        # File clinvar_vcf_uri
        # File omim_uri
        # String mpc_ht_uri
        # String loeuf_v2_uri
        # String loeuf_v4_uri

        String cohort_prefix
        String hail_docker
        String vep_hail_docker
        String sv_base_mini_docker
        
        String vep_annotate_hail_python_script
        String split_vcf_hail_script

        String genome_build='GRCh38'
        Boolean split_by_chromosome
        Boolean split_into_shards 
        Boolean merge_split_vcf
        Boolean reannotate_ac_af=false
        Int shards_per_chunk=10  # combine pre-sharded VCFs
        
        Array[File]? vcf_shards  # if scatterVCF.wdl already run before VEP
        
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_vep_annotate
        RuntimeAttr? runtime_attr_annotate_add_genotypes
    }

    if (defined(vcf_shards)) {
        String file_ = select_first([select_first([vcf_shards])[0]])
    }
    String file = select_first([file_, vcf_file])

    # input is not a single VCF file, so merge shards in chunks, then run VEP on merged chunks
    if (merge_split_vcf) { 
        # combine pre-sharded VCFs into chunks
        call mergeSplitVCF.splitFile as splitFile {
            input:
                file=file,
                shards_per_chunk=shards_per_chunk,
                cohort_prefix=cohort_prefix,
                hail_docker=vep_hail_docker
        }
        scatter (chunk_file in splitFile.chunks) {
            call mergeVCFs.mergeVCFs as mergeVCFs {
                input:
                    vcf_files=read_lines(chunk_file),
                    sv_base_mini_docker=sv_base_mini_docker,
                    cohort_prefix=basename(chunk_file),
                    runtime_attr_override=runtime_attr_merge_vcfs
            }
            
            call RemoveRawMutectCalls as RemoveRawMutectCallsMergedShards {
                input:
                vcf=mergeVCFs.merged_vcf_file,
                sv_base_mini_docker=sv_base_mini_docker
            }

            call vepAnnotate as vepAnnotateMergedShards {
                input:
                    vcf_file=RemoveRawMutectCallsMergedShards.out,
                    vep_annotate_hail_python_script=vep_annotate_hail_python_script,
                    top_level_fa=top_level_fa,
                    # human_ancestor_fa=human_ancestor_fa,
                    # human_ancestor_fa_fai=human_ancestor_fa_fai,
                    # gerp_conservation_scores=gerp_conservation_scores,
                    ref_vep_cache=ref_vep_cache,
                    alpha_missense_file=alpha_missense_file,
                    alpha_missense_file_idx=alpha_missense_file+'.tbi',
                    eve_data=eve_data,
                    eve_data_idx=eve_data+'.tbi',
                    vep_hail_docker=vep_hail_docker,
                    reannotate_ac_af=reannotate_ac_af,
                    genome_build=genome_build,
                    runtime_attr_override=runtime_attr_vep_annotate
            }

            call addHGVScGeneSymbolField as addHGVScGeneSymbolFieldMergedShards {
                input:
                vcf_file=vepAnnotateMergedShards.vep_vcf_file,
                genome_build=genome_build,
                vep_hail_docker=vep_hail_docker,
                runtime_attr_override=runtime_attr_vep_annotate
            }

            call annotateExtra.addGenotypes as addGenotypesMergedShards {
                input:
                annot_vcf_file=addHGVScGeneSymbolFieldMergedShards.output_vcf,
                vcf_file=mergeVCFs.merged_vcf_file,
                hail_docker=hail_docker,
                genome_build=genome_build,
                runtime_attr_override=runtime_attr_annotate_add_genotypes
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
                    genome_build=genome_build,
                    hail_docker=hail_docker,
                    sv_base_mini_docker=sv_base_mini_docker,
                    split_by_chromosome=split_by_chromosome,
                    split_into_shards=split_into_shards,
            }
        }
        Array[File] vcf_shards_ = select_first([scatterVCF.vcf_shards, vcf_shards])
    
        scatter (vcf_shard in vcf_shards_) {
            call RemoveRawMutectCalls {
                input:
                vcf=vcf_shard,
                sv_base_mini_docker=sv_base_mini_docker
            }

            call vepAnnotate {
                input:
                    vcf_file=RemoveRawMutectCalls.out,
                    vep_annotate_hail_python_script=vep_annotate_hail_python_script,
                    top_level_fa=top_level_fa,
                    # human_ancestor_fa=human_ancestor_fa,
                    # human_ancestor_fa_fai=human_ancestor_fa_fai,
                    # gerp_conservation_scores=gerp_conservation_scores,
                    ref_vep_cache=ref_vep_cache,
                    alpha_missense_file=alpha_missense_file,
                    alpha_missense_file_idx=alpha_missense_file+'.tbi',
                    eve_data=eve_data,
                    eve_data_idx=eve_data+'.tbi',
                    vep_hail_docker=vep_hail_docker,
                    reannotate_ac_af=reannotate_ac_af,
                    genome_build=genome_build,
                    runtime_attr_override=runtime_attr_vep_annotate
            }

            call addHGVScGeneSymbolField {
                input:
                    vcf_file=vepAnnotate.vep_vcf_file,
                    genome_build=genome_build,
                    vep_hail_docker=vep_hail_docker,
                    runtime_attr_override=runtime_attr_vep_annotate
            }

            call annotateExtra.addGenotypes as addGenotypes {
                input:
                annot_vcf_file=addHGVScGeneSymbolField.output_vcf,
                vcf_file=vcf_shard,
                hail_docker=hail_docker,
                genome_build=genome_build,
                runtime_attr_override=runtime_attr_annotate_add_genotypes
            }

        }
    }

    Array[File] vep_vcf_files_ = select_first([addGenotypesMergedShards.combined_vcf_file, addGenotypes.combined_vcf_file])
    Array[File] vep_vcf_idx_ = select_first([addGenotypesMergedShards.combined_vcf_idx, addGenotypes.combined_vcf_idx])

    output {
        Array[File] vep_vcf_files = vep_vcf_files_
        Array[File] vep_vcf_idx = vep_vcf_idx_
    }
}   

task RemoveRawMutectCalls {
    input {
        File vcf
        String sv_base_mini_docker
    }

    RuntimeAttr runtime_default = object {
                                      mem_gb: 3,
                                      disk_gb: 10,
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }

    RuntimeAttr runtime_override = runtime_default
    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])

    String prefix = basename(vcf, ".vcf.gz")

    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -eou pipefail
        # mutect calls are the first sample
        bcftools query -l ~{vcf} | awk 'NR > 1' > samples.txt
        bcftools view ~{vcf} -S samples.txt -Oz -o ~{prefix}.fm.vcf.gz

    >>>

    output {
        File out = "~{prefix}.fm.vcf.gz"
    }
}

task vepAnnotate {
    input {
        File vcf_file
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
        String vep_annotate_hail_python_script
        Boolean reannotate_ac_af
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB") + size(ref_vep_cache, "GB")# + size(gerp_conservation_scores, "GB")
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

        curl ~{vep_annotate_hail_python_script} > vep_annotate.py
        proj_id=$(gcloud config get-value project)
        python3.9 vep_annotate.py -i ~{vcf_file} -o ~{vep_annotated_vcf_name} --cores ~{cpu_cores} --mem ~{memory} \
        --reannotate-ac-af ~{reannotate_ac_af} --build ~{genome_build} --project-id $proj_id
        cp $(ls . | grep hail*.log) hail_log.txt
        bcftools index -t ~{vep_annotated_vcf_name}
    >>>

    output {
        File vep_vcf_file = vep_annotated_vcf_name
        File vep_vcf_idx = vep_annotated_vcf_name + '.tbi'
        File hail_log = "hail_log.txt"
    }
}

task addHGVScGeneSymbolField {
    input {
        File vcf_file
        String genome_build
        String vep_hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
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

    String output_filename = basename(vcf_file, '.vcf.bgz') + '_HGVSc_symbol.vcf.bgz'

    command <<<
        set -euo pipefail

        cat <<EOF > add_hgvsc.py
        import hail as hl
        import os
        import sys

        vcf_file = sys.argv[1]
        build = sys.argv[2]
        output_filename = sys.argv[3]

        hl.init()

        mt = hl.import_vcf(vcf_file, reference_genome=build, force_bgz=True)

        header = hl.get_vcf_metadata(vcf_file)
        csq_columns = header['info']['CSQ']['Description'].split('Format: ')[1].split('|')

        mt = mt.annotate_rows(vep=mt.info)
        transcript_consequences = mt.vep.CSQ.map(lambda x: x.split('\|'))

        transcript_consequences_strs = transcript_consequences.map(lambda x: hl.if_else(hl.len(x)>1, hl.struct(**
                                                               {col: x[i] if col!='Consequence' else x[i].split('&')  
                                                                for i, col in enumerate(csq_columns)}), 
                                                                hl.struct(**{col: hl.missing('str') if col!='Consequence' else hl.array([hl.missing('str')])  
                                                                for i, col in enumerate(csq_columns)})))

        mt = mt.annotate_rows(vep=mt.vep.annotate(transcript_consequences=transcript_consequences_strs))
        mt = mt.annotate_rows(vep=mt.vep.select('transcript_consequences'))

        mt_by_gene = mt.explode_rows(mt.vep.transcript_consequences)

        mt_by_gene = mt_by_gene.annotate_rows(vep=mt_by_gene.vep.annotate(
            transcript_consequences=mt_by_gene.vep.transcript_consequences.annotate(
                HGVSc_symbol= hl.str(', ').join(
                            hl.array([mt_by_gene.vep.transcript_consequences.SYMBOL, 
                                     hl.if_else(mt_by_gene.vep.transcript_consequences.HGVSc!='', mt_by_gene.vep.transcript_consequences.HGVSc.split(':')[1], hl.missing('str')),
                                     hl.if_else(mt_by_gene.vep.transcript_consequences.HGVSp!='', mt_by_gene.vep.transcript_consequences.HGVSp.split(':')[1], hl.missing('str'))])
            .filter(lambda x: hl.is_defined(x))
                                        ))))

        csq_fields_str = '|'.join(csq_columns + ['HGVSc_symbol'])

        mt_by_gene = (mt_by_gene.group_rows_by(mt_by_gene.locus, mt_by_gene.alleles)
            .aggregate_rows(vep = hl.agg.collect(mt_by_gene.vep))).result()

        fields = list(mt_by_gene.vep.transcript_consequences[0])
        new_csq = mt_by_gene.vep.transcript_consequences.scan(lambda i, j: 
                                              hl.str('|').join(hl.array([i]))
                                              +','+hl.str('|').join(hl.array([j[col] if col!='Consequence' else 
                                                                          hl.str('&').join(j[col]) 
                                                                          for col in list(fields)])), '')[-1][1:]
        mt_by_gene = mt_by_gene.annotate_rows(CSQ=new_csq)

        mt = mt.annotate_rows(info=mt.info.annotate(CSQ=mt_by_gene.rows()[mt.row_key].CSQ))
        mt = mt.drop('vep')

        header['info']['CSQ'] = {'Description': csq_fields_str, 'Number': '.', 'Type': 'String'}

        hl.export_vcf(dataset=mt, output=output_filename, metadata=header)
        EOF

        python3.9 add_hgvsc.py ~{vcf_file} ~{genome_build} ~{output_filename}

        /opt/vep/bcftools/bcftools index -t ~{output_filename}
    >>>

    output {
        File output_vcf = output_filename
        File output_vcf_idx = output_filename + '.tbi'
    }
}