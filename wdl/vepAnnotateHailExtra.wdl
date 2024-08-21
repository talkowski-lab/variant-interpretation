version 1.0
    
import "scatterVCF.wdl" as scatterVCF
import "mergeSplitVCF.wdl" as mergeSplitVCF
import "mergeVCFs.wdl" as mergeVCFs
import "annotateNonCoding.wdl" as annotateNonCoding

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow vepAnnotateHailExtra {

    input {
        Array[File] vep_vcf_files

        File revel_file
        File clinvar_vcf_uri
        File omim_uri
        String mpc_ht_uri
        String loeuf_v2_uri
        String loeuf_v4_uri

        String cohort_prefix
        String hail_docker
        String sv_base_mini_docker
        
        String vep_annotate_hail_extra_python_script
        String split_vcf_hail_script

        String genome_build='GRCh38'

        String spliceAI_snv_uri='NA'
        String spliceAI_indel_uri='NA'       
        String noncoding_bed='NA'
        String gene_list='NA'

        RuntimeAttr? runtime_attr_annotate_noncoding      
        RuntimeAttr? runtime_attr_annotate_extra
        RuntimeAttr? runtime_attr_annotate_spliceAI
    }

    scatter (vcf_shard in vep_vcf_files) {

        if (noncoding_bed!='NA') {
            call annotateNonCoding.annotateFromBed as annotateNonCoding {
                input:
                vcf_file=vcf_shard,
                noncoding_bed=select_first([noncoding_bed]),
                hail_docker=hail_docker,
                runtime_attr_override=runtime_attr_annotate_noncoding
            }
        }
        call annotateExtra {
            input:
                vcf_file=select_first([annotateNonCoding.noncoding_vcf, vcf_shard]),
                vep_annotate_hail_extra_python_script=vep_annotate_hail_extra_python_script,
                loeuf_v2_uri=loeuf_v2_uri,
                loeuf_v4_uri=loeuf_v4_uri,
                revel_file=revel_file,
                revel_file_idx=revel_file+'.tbi',
                clinvar_vcf_uri=clinvar_vcf_uri,
                omim_uri=omim_uri,
                gene_list=select_first([gene_list, 'NA']),
                mpc_ht_uri=mpc_ht_uri,
                hail_docker=hail_docker,
                genome_build=genome_build,
                runtime_attr_override=runtime_attr_annotate_extra
        }

        if ((spliceAI_snv_uri!='NA') && (spliceAI_indel_uri!='NA')) {
            call annotateSpliceAI {
                input:
                vcf_file=annotateExtra.annot_vcf_file,
                spliceAI_snv_uri=spliceAI_snv_uri,
                spliceAI_indel_uri=spliceAI_indel_uri,
                genome_build=genome_build,
                hail_docker=hail_docker,
                runtime_attr_override=runtime_attr_annotate_spliceAI
            }
        }
        File annot_vcf_file = select_first([annotateSpliceAI.annot_vcf_file, annotateExtra.annot_vcf_file])
        File annot_vcf_idx_ = select_first([annotateSpliceAI.annot_vcf_idx, annotateExtra.annot_vcf_idx])
    }

    output {
        Array[File] annot_vcf_files = annot_vcf_file
        Array[File] annot_vcf_idx = annot_vcf_idx_
    }
}   

task annotateExtra {
    input {
        File vcf_file
        String loeuf_v2_uri
        String loeuf_v4_uri
        File revel_file
        File revel_file_idx
        File clinvar_vcf_uri
        File omim_uri
        
        String gene_list
        String mpc_ht_uri

        String hail_docker
        String genome_build
        String vep_annotate_hail_extra_python_script
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
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String filename = basename(vcf_file)
    String prefix = if (sub(filename, "\\.gz", "")!=filename) then basename(vcf_file, ".vcf.gz") else basename(vcf_file, ".vcf.bgz")
    String vep_annotated_vcf_name = "~{prefix}.annot.vcf.bgz"

    command <<<
        curl ~{vep_annotate_hail_extra_python_script} > annotate.py
        python3 annotate.py -i ~{vcf_file} -o ~{vep_annotated_vcf_name} --cores ~{cpu_cores} --mem ~{memory} \
        --build ~{genome_build} --loeuf-v2 ~{loeuf_v2_uri} --loeuf-v4 ~{loeuf_v4_uri} \
        --mpc ~{mpc_ht_uri} --clinvar ~{clinvar_vcf_uri} --omim ~{omim_uri} \
        --revel ~{revel_file} --genes ~{gene_list} 
        cp $(ls . | grep hail*.log) hail_log.txt
    >>>

    output {
        File annot_vcf_file = vep_annotated_vcf_name
        File annot_vcf_idx = vep_annotated_vcf_name + '.tbi'
        File hail_log = "hail_log.txt"
    }
}

task annotateSpliceAI {
    input {
        File vcf_file
        String spliceAI_snv_uri
        String spliceAI_indel_uri

        String hail_docker
        String genome_build
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
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String filename = basename(vcf_file)
    String prefix = if (sub(filename, "\\.gz", "")!=filename) then basename(vcf_file, ".vcf.gz") else basename(vcf_file, ".vcf.bgz")
    String vep_annotated_vcf_name = "~{prefix}.SpliceAI.annot.vcf.bgz"


    command <<<
    cat <<EOF > annotate.py
    from pyspark.sql import SparkSession
    import hail as hl
    import numpy as np
    import pandas as pd
    import sys
    import ast
    import os
    import json
    import argparse

    parser = argparse.ArgumentParser(description='Parse arguments')
    parser.add_argument('-i', dest='vcf_file', help='Input VCF file')
    parser.add_argument('-o', dest='vep_annotated_vcf_name', help='Output filename')
    parser.add_argument('--cores', dest='cores', help='CPU cores')
    parser.add_argument('--mem', dest='mem', help='Memory')
    parser.add_argument('--build', dest='build', help='Genome build')
    parser.add_argument('--spliceAI-snv', dest='spliceAI_snv_uri', help='SpliceAI scores SNV HT')
    parser.add_argument('--spliceAI-indel', dest='spliceAI_indel_uri', help='SpliceAI scores Indel HT')

    args = parser.parse_args()

    vcf_file = args.vcf_file
    vep_annotated_vcf_name = args.vep_annotated_vcf_name
    cores = args.cores  # string
    mem = int(np.floor(float(args.mem)))
    build = args.build
    spliceAI_snv_uri = args.spliceAI_snv_uri
    spliceAI_indel_uri = args.spliceAI_indel_uri

    hl.init(min_block_size=128, 
            local=f"local[{cores}]", 
            spark_conf={
                # "spark.executor.cores": '1', 
                #         "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                #         "spark.driver.cores": cores,
                        "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                        "spark.driver.memoryOverheadFactor": '0.6',
                        # "spark.executor.memoryOverheadFactor": '0.6',
            #             'spark.hadoop.fs.gs.requester.pays.mode': 'CUSTOM',
            #             'spark.hadoop.fs.gs.requester.pays.buckets': 'hail-datasets-us-central1',
            #             'spark.hadoop.fs.gs.requester.pays.project.id': gcp_project,
                        }, 
            tmp_dir="tmp", local_tmpdir="tmp",
                        )

    header = hl.get_vcf_metadata(vcf_file) 
    mt = hl.import_vcf(vcf_file, force_bgz=True, array_elements_required=False, call_fields=[], reference_genome=build)

    csq_columns = header['info']['CSQ']['Description'].split('Format: ')[1].split('|')
    # split VEP CSQ string
    mt = mt.annotate_rows(vep=mt.info)
    transcript_consequences = mt.vep.CSQ.map(lambda x: x.split('\|'))

    transcript_consequences_strs = transcript_consequences.map(lambda x: hl.if_else(hl.len(x)>1, hl.struct(**
                                                        {col: x[i] if col!='Consequence' else x[i].split('&')  
                                                            for i, col in enumerate(csq_columns)}), 
                                                            hl.struct(**{col: hl.missing('str') if col!='Consequence' else hl.array([hl.missing('str')])  
                                                            for i, col in enumerate(csq_columns)})))

    mt = mt.annotate_rows(vep=mt.vep.annotate(transcript_consequences=transcript_consequences_strs))
    mt = mt.annotate_rows(vep=mt.vep.select('transcript_consequences'))

    # annotate SpliceAI scores
    mt_by_transcript = mt.explode_rows(mt.vep.transcript_consequences)
    mt_by_locus_and_gene = mt_by_transcript.key_rows_by('locus', 'alleles', mt_by_transcript.vep.transcript_consequences.SYMBOL)

    snv_ht = hl.read_table(spliceAI_snv_uri)
    indel_ht = hl.read_table(spliceAI_indel_uri)
    spliceAI_ht = snv_ht.union(indel_ht)
    # leave out ALLELE/SYMBOL because redundant
    fields = 'ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL'.split('|')[2:]  
    mt_by_locus_and_gene = mt_by_locus_and_gene.annotate_rows(SpliceAI_raw=spliceAI_ht[mt_by_locus_and_gene.row_key].SpliceAI)
    mt_by_locus_and_gene = mt_by_locus_and_gene.annotate_rows(vep=mt_by_locus_and_gene.vep.annotate(
        transcript_consequences=(mt_by_locus_and_gene.vep.transcript_consequences.annotate(
            **{field: hl.if_else(hl.is_defined(mt_by_locus_and_gene.SpliceAI_raw), 
                                mt_by_locus_and_gene.SpliceAI_raw.split('=')[1].split('\|')[i+2], '') 
            for i, field in enumerate(fields)}))))
    csq_fields_str = 'Format: ' + header['info']['CSQ']['Description'].split('Format: ')[1] + '|'.join([''] + fields)

    mt_by_gene = mt_by_locus_and_gene
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
    header['info']['REVEL'] = {'Description': 'REVEL scores.', 'Number': '.', 'Type': 'String'}

    hl.export_vcf(dataset=mt, output=vep_annotated_vcf_name, metadata=header, tabix=True)
    EOF
    python3 annotate.py -i ~{vcf_file} -o ~{vep_annotated_vcf_name} --cores ~{cpu_cores} --mem ~{memory} \
    --build ~{genome_build} --spliceAI-snv ~{spliceAI_snv_uri} --spliceAI-indel ~{spliceAI_indel_uri}
    cp $(ls . | grep hail*.log) hail_log.txt
    >>>

    output {
        File annot_vcf_file = vep_annotated_vcf_name
        File annot_vcf_idx = vep_annotated_vcf_name + '.tbi'
        File hail_log = "hail_log.txt"
    }
}