from pyspark.sql import SparkSession
import hail as hl
import numpy as np
import pandas as pd
import sys
import ast
import os
import json
import argparse
import datetime

parser = argparse.ArgumentParser(description='Parse arguments')
parser.add_argument('-i', dest='ht_uri', help='Input HT')
parser.add_argument('--bucket-id', dest='bucket_id', help='Bucket ID')
parser.add_argument('--cores', dest='cores', help='CPU cores')
parser.add_argument('--mem', dest='mem', help='Memory')
parser.add_argument('--mpc', dest='mpc_ht_uri', help='MPC scores HT')
parser.add_argument('--clinvar', dest='clinvar_vcf_uri', help='ClinVar VCF')
parser.add_argument('--omim', dest='omim_uri', help='OMIM file')
parser.add_argument('--revel', dest='revel_file', help='REVEL file')
parser.add_argument('--loeuf-v2', dest='loeuf_v2_uri', help='LOEUF scores from gnomAD v2.1.1')
parser.add_argument('--loeuf-v4', dest='loeuf_v4_uri', help='LOEUF scores from gnomAD v4.1')
parser.add_argument('--genes', dest='gene_list', help='OPTIONAL: Gene list txt file')
parser.add_argument('--build', dest='build', help='Genome build')
parser.add_argument('--project-id', dest='project_id', help='Google Project ID')

args = parser.parse_args()

ht_uri = args.ht_uri
bucket_id = args.bucket_id
cores = args.cores  # string
mem = int(np.floor(float(args.mem)))
build = args.build
gcp_project = args.project_id
mpc_ht_uri = args.mpc_ht_uri
clinvar_vcf_uri = args.clinvar_vcf_uri
omim_uri = args.omim_uri
revel_file = args.revel_file
loeuf_v2_uri = args.loeuf_v2_uri
loeuf_v4_uri = args.loeuf_v4_uri
gene_list = args.gene_list
gcp_project = args.project_id

hl.init(min_block_size=128, 
        spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{int(np.floor(mem*0.4))}g",
        #             'spark.hadoop.fs.gs.requester.pays.mode': 'CUSTOM',
        #             'spark.hadoop.fs.gs.requester.pays.buckets': 'hail-datasets-us-central1',
        #             'spark.hadoop.fs.gs.requester.pays.project.id': gcp_project,
                    }, 
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

ht = hl.read_table(ht_uri)

# run VEP
ht = hl.vep(ht, config='vep_config.json', csq=True, tolerate_parse_error=True)

## TODO: set up after VEP finished
# FROM EXTRA
# annotate MPC
mpc = hl.read_table(mpc_ht_uri).key_by('locus','alleles')
ht = ht.annotate(MPC=mpc[ht.locus, ht.alleles].mpc)
        
# annotate ClinVar
if build=='GRCh38':
    clinvar_vcf = hl.import_vcf(clinvar_vcf_uri,
                            reference_genome='GRCh38',
                            force_bgz=clinvar_vcf_uri.split('.')[-1] in ['gz', 'bgz'])
    ht = ht.annotate(CLNSIG=clinvar_vcf.rows()[ht.key].info.CLNSIG,
                                                  CLNREVSTAT=clinvar_vcf.rows()[ht.key].info.CLNREVSTAT)

# annotate REVEL
revel_ht = hl.import_table(revel_file, force_bgz=True)
revel_ht = revel_ht.annotate(chr='chr'+revel_ht['#chr']) 
build_chr = 'chr' if build=='GRCh38' else '#chr'
build_pos = 'grch38_pos' if build=='GRCh38' else 'hg19_pos'
revel_ht = revel_ht.annotate(locus=hl.locus(revel_ht[build_chr], hl.int(revel_ht[build_pos]), build),
                 alleles=hl.array([revel_ht.ref, revel_ht.alt]))
revel_ht = revel_ht.key_by('locus', 'alleles')
ht = ht.annotate(REVEL=revel_ht[ht.key].REVEL)

# split VEP CSQ string
transcript_consequences = ht.vep.map(lambda x: x.split('\|'))

csq_columns = hl.eval(ht.vep_csq_header).split('Format: ')[1].split('|')
transcript_consequences_strs = transcript_consequences.map(lambda x: hl.if_else(hl.len(x)>1, hl.struct(**
                                                       {col: x[i] if col!='Consequence' else x[i].split('&')  
                                                        for i, col in enumerate(csq_columns)}), 
                                                        hl.struct(**{col: hl.missing('str') if col!='Consequence' else hl.array([hl.missing('str')])  
                                                        for i, col in enumerate(csq_columns)})))

ht = ht.annotate(vep=hl.Struct(**{'transcript_consequences': transcript_consequences_strs}))

# annotate LOEUF from gnomAD
loeuf_v2_ht = hl.read_table(loeuf_v2_uri).key_by('transcript')
loeuf_v4_ht = hl.read_table(loeuf_v4_uri).key_by('transcript')
ht_by_transcript = ht.explode(ht.vep.transcript_consequences)
ht_by_transcript = ht_by_transcript.key_by(ht_by_transcript.vep.transcript_consequences.Feature)
ht_by_transcript = ht_by_transcript.annotate(vep=ht_by_transcript.vep.annotate(
    transcript_consequences=ht_by_transcript.vep.transcript_consequences.annotate(
        LOEUF_v2=hl.if_else(hl.is_defined(loeuf_v2_ht[ht_by_transcript.key]), loeuf_v2_ht[ht_by_transcript.key]['oe_lof_upper'], ''),
        LOEUF_v2_decile=hl.if_else(hl.is_defined(loeuf_v2_ht[ht_by_transcript.key]), loeuf_v2_ht[ht_by_transcript.key]['oe_lof_upper_bin'], ''),
        LOEUF_v4=hl.if_else(hl.is_defined(loeuf_v4_ht[ht_by_transcript.key]), loeuf_v4_ht[ht_by_transcript.key]['lof.oe_ci.upper'], ''),
        LOEUF_v4_decile=hl.if_else(hl.is_defined(loeuf_v4_ht[ht_by_transcript.key]), loeuf_v4_ht[ht_by_transcript.key]['lof.oe_ci.upper_bin_decile'], '')
        )
    )
)
csq_fields_str = hl.eval(ht.vep_csq_header) + '|'.join(['', 'LOEUF_v2', 'LOEUF_v2_decile', 'LOEUF_v4', 'LOEUF_v4_decile'])

# annotate OMIM
omim = hl.import_table(omim_uri).key_by('approvedGeneSymbol')
ht_by_gene = ht_by_transcript.key_by(ht_by_transcript.vep.transcript_consequences.SYMBOL)
ht_by_gene = ht_by_gene.annotate(vep=ht_by_gene.vep.annotate(
    transcript_consequences=ht_by_gene.vep.transcript_consequences.annotate(    
    OMIM_MIM_number=hl.if_else(hl.is_defined(omim[ht_by_gene.key]), omim[ht_by_gene.key].mimNumber, ''),
    OMIM_inheritance_code=hl.if_else(hl.is_defined(omim[ht_by_gene.key]), omim[ht_by_gene.key].inheritance_code, '')))
)
csq_fields_str = csq_fields_str + '|'.join([''] + ['OMIM_MIM_number', 'OMIM_inheritance_code'])

# OPTIONAL: annotate with gene list, if provided
if gene_list!='NA':
    genes = pd.read_csv(gene_list, sep='\t', header=None)[0].tolist()
    gene_list_name = os.path.basename(gene_list)
    ht_by_gene = ht_by_gene.annotate(vep=ht_by_gene.vep.annotate(
    transcript_consequences=ht_by_gene.vep.transcript_consequences.annotate(    
        gene_list=hl.if_else(hl.array(genes).contains(ht_by_gene.key.SYMBOL), gene_list_name, '')))
    )
    csq_fields_str = csq_fields_str + '|gene_list'

ht_by_gene = (ht_by_gene.group_by(ht_by_gene.locus, ht_by_gene.alleles)
    .aggregate(transcript_consequences = hl.agg.collect(ht_by_gene.vep.transcript_consequences)))

ht = ht.annotate(vep=hl.Struct(**{'transcript_consequences': ht_by_gene[ht.key].transcript_consequences}))
ht = ht.drop('vep_csq_header')
ht = ht.annotate(vep_csq_header=csq_fields_str)

prefix = os.path.basename(ht_uri).split('.ht')[0]
filename = f"{bucket_id}/hail/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))}/{prefix}.vep.ht"
pd.Series([filename]).to_csv('ht_uri.txt', index=False, header=None)
ht.write(filename)