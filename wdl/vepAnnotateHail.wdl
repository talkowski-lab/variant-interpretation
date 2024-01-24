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
        File file
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
        String vep_hail_docker
        String sv_base_mini_docker
        Boolean split_by_chromosome
        Boolean split_into_shards 
        Int compression_level=3
        Int shards_per_chunk=10  # combine pre-sharded VCFs
        Int n_shards=0 
        Array[File]? vcf_shards  # if scatter-vcf run before VEP
        Int? thread_num_override
        Float? input_disk_scale_merge_chunk
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_split_vcf
        RuntimeAttr? runtime_attr_vep_annotate
    }

    String filename = basename(file)
    if (sub(filename, ".vcf.gz", "")==filename) {  # input is not a single VCF file
        Boolean merge_shards=true
        # combine pre-sharded VCFs into chunks
        call splitFile {
            input:
                file=file,
                shards_per_chunk=shards_per_chunk,
                cohort_prefix=cohort_prefix,
                vep_hail_docker=vep_hail_docker
        }
    }

    # shard the VCF (if not already sharded)
    if (split_by_chromosome) {
        call splitByChromosome {
            input:
                vcf_file=file,
                sv_base_mini_docker=sv_base_mini_docker,
                thread_num_override=thread_num_override,
                compression_level=compression_level,
                runtime_attr_override=runtime_attr_split_vcf
        }
    }
    if (split_into_shards) {
        # if already split into chromosomes, shard further
        if (defined(splitByChromosome.shards)) {
            scatter (chrom_shard in select_first([splitByChromosome.shards])) {
                call scatterVCF as scatterChromosomes {
                    input:
                        vcf_file=chrom_shard,
                        split_vcf_hail_script=split_vcf_hail_script,
                        n_shards=select_first([n_shards]),
                        vep_hail_docker=vep_hail_docker,
                        thread_num_override=thread_num_override,
                        compression_level=compression_level,
                        runtime_attr_override=runtime_attr_split_vcf
                }
            }
            Array[File] chromosome_shards = flatten(scatterChromosomes.shards)
        }

        if (!defined(splitByChromosome.shards)) {
            call scatterVCF {
                input:
                    vcf_file=file,
                    split_vcf_hail_script=split_vcf_hail_script,
                    n_shards=select_first([n_shards]),
                    vep_hail_docker=vep_hail_docker,
                    thread_num_override=thread_num_override,
                    compression_level=compression_level,
                    runtime_attr_override=runtime_attr_split_vcf
            }
        }
    }
    
    Array[File] vcf_shards_ = select_first([scatterVCF.shards, chromosome_shards, 
                                        splitByChromosome.shards, splitFile.chunks, vcf_shards, [file]])

    # if split into chunks, merge shards in chunks, then run VEP on merged chunks
    if (defined(merge_shards)) {
        scatter (chunk_file in vcf_shards_) {
            call mergeVCFs {
                input:
                    vcf_files=read_lines(chunk_file),
                    sv_base_mini_docker=sv_base_mini_docker,
                    cohort_prefix=basename(chunk_file),
                    merge_or_concat='concat',
                    input_disk_scale=input_disk_scale_merge_chunk,
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
       
    if (!defined(merge_shards)) {
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

task mergeVCFs {
    input {
        Array[File] vcf_files
        String sv_base_mini_docker
        String cohort_prefix
        String merge_or_concat 
        Float? input_disk_scale
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_files, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale_ = select_first([input_disk_scale, 5.0])
    
    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale_),
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

    String merged_vcf_name="~{cohort_prefix}.merged.vcf.gz"
    String sorted_vcf_name="~{cohort_prefix}.merged.sorted.vcf.gz"

    String merge_or_concat_new = if merge_or_concat == 'concat' then 'concat -n'  else merge_or_concat
    command <<<
        set -euo pipefail
        VCFS="~{write_lines(vcf_files)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        bcftools ~{merge_or_concat_new} --no-version -Oz --file-list vcfs_sorted.list --output ~{merged_vcf_name}
        bcftools sort ~{merged_vcf_name} --output ~{sorted_vcf_name}
        bcftools index -t ~{sorted_vcf_name}
    >>>

    output {
        File merged_vcf_file=sorted_vcf_name
        File merged_vcf_idx=sorted_vcf_name + ".tbi"
    }
}

task splitFile {
    input {
        File file
        Int shards_per_chunk
        String cohort_prefix
        String vep_hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(file, "GB")
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
        docker: vep_hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        split -l ~{shards_per_chunk} ~{file} -a 4 -d "~{cohort_prefix}.shard."
    >>>

    output {
        Array[File] chunks = glob("~{cohort_prefix}.*")
    }
}

task splitByChromosome { 
    input {
        File vcf_file
        String sv_base_mini_docker
        Int? thread_num_override
        Int? compression_level
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
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
    Int thread_num = select_first([thread_num_override,runtime_override.cpu_cores])
    String compression_str = if defined(compression_level) then "-Oz~{compression_level}" else "-Oz"

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
        awk '$0="~{prefix}."$0' chr_list.txt > filenames.txt
        paste chr_list.txt filenames.txt > chr_filenames.txt
        bcftools +scatter ~{vcf_file} ~{compression_str} -o . --threads ~{thread_num} -S chr_filenames.txt 
    >>>

    output {
        Array[File] shards = glob("~{prefix}.chr*.vcf.gz")
        Array[String] shards_string = glob("~{prefix}.chr*.vcf.gz")
    }
}

task scatterVCF {
    input {
        File vcf_file
        String split_vcf_hail_script
        Int n_shards
        String vep_hail_docker
        Int? compression_level
        Int? thread_num_override
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
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
        curl ~{split_vcf_hail_script} > split_vcf.py
        python3.9 split_vcf.py ~{vcf_file} ~{n_shards} ~{prefix} ~{cpu_cores} ~{memory}
        for file in $(ls ~{prefix}.vcf.bgz | grep '.bgz'); do
            shard_num=$(echo $file | cut -d '-' -f2);
            mv ~{prefix}.vcf.bgz/$file ~{prefix}.shard_"$shard_num".vcf.bgz
        done
    >>>
    output {
        Array[File] shards = glob("~{prefix}.shard_*.vcf.bgz")
        Array[String] shards_string = glob("~{prefix}.shard_*.vcf.bgz")
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
    String prefix = if (sub(vcf_file, ".gz", "")!=filename) then basename(vcf_file, ".vcf.gz") else basename(vcf_file, ".vcf.bgz")
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