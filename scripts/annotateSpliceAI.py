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