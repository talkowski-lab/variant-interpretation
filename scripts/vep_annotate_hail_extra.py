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
parser.add_argument('--mpc', dest='mpc_ht_uri', help='MPC scores HT')
parser.add_argument('--clinvar', dest='clinvar_vcf_uri', help='ClinVar VCF')
parser.add_argument('--omim', dest='omim_uri', help='OMIM file')
parser.add_argument('--revel', dest='revel_file', help='REVEL file')
parser.add_argument('--loeuf-v2', dest='loeuf_v2_uri', help='LOEUF scores from gnomAD v2.1.1')
parser.add_argument('--loeuf-v4', dest='loeuf_v4_uri', help='LOEUF scores from gnomAD v4.1')
parser.add_argument('--spliceAI-snv', dest='spliceAI_snv_uri', help='SpliceAI scores SNV HT')
parser.add_argument('--spliceAI-indel', dest='spliceAI_indel_uri', help='SpliceAI scores Indel HT')
parser.add_argument('--genes', dest='gene_list', help='OPTIONAL: Gene list txt file')
parser.add_argument('--noncoding-bed', dest='noncoding_bed', help='OPTIONAL: non-coding annotations bed file')
parser.add_argument('--project-id', dest='project_id', help='Google Project ID')

args = parser.parse_args()

vcf_file = args.vcf_file
vep_annotated_vcf_name = args.vep_annotated_vcf_name
cores = args.cores  # string
mem = int(np.floor(float(args.mem)))
build = args.build
mpc_ht_uri = args.mpc_ht_uri
clinvar_vcf_uri = args.clinvar_vcf_uri
omim_uri = args.omim_uri
revel_file = args.revel_file
loeuf_v2_uri = args.loeuf_v2_uri
loeuf_v4_uri = args.loeuf_v4_uri
spliceAI_snv_uri = args.spliceAI_snv_uri
spliceAI_indel_uri = args.spliceAI_indel_uri
gene_list = args.gene_list
noncoding_bed = args.noncoding_bed
gcp_project = args.project_id

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

# OPTIONAL: annotate non-coding 
if noncoding_bed!='NA':
    bed = hl.import_bed(noncoding_bed, reference_genome=build, skip_invalid_intervals=True)
    mt = mt.annotate_rows(info=mt.info.annotate(PREDICTED_NONCODING=bed[mt.locus].target))
    mt.checkpoint('noncoding_annot.mt')

# annotate MPC
mpc = hl.read_table(mpc_ht_uri).key_by('locus','alleles')
mt = mt.annotate_rows(info = mt.info.annotate(MPC=mpc[mt.locus, mt.alleles].mpc))
        
# annotate ClinVar
if build=='GRCh38':
    clinvar_vcf = hl.import_vcf(clinvar_vcf_uri,
                            reference_genome='GRCh38',
                            force_bgz=clinvar_vcf_uri.split('.')[-1] in ['gz', 'bgz'])
    mt = mt.annotate_rows(info = mt.info.annotate(CLNSIG=clinvar_vcf.rows()[mt.row_key].info.CLNSIG,
                                                  CLNREVSTAT=clinvar_vcf.rows()[mt.row_key].info.CLNREVSTAT))

# annotate REVEL
revel_ht = hl.import_table(revel_file, force_bgz=True)
revel_ht = revel_ht.annotate(chr='chr'+revel_ht['#chr']) 
build_chr = 'chr' if build=='GRCh38' else '#chr'
build_pos = 'grch38_pos' if build=='GRCh38' else 'hg19_pos'
revel_ht = revel_ht.annotate(locus=hl.locus(revel_ht[build_chr], hl.int(revel_ht[build_pos]), build),
                 alleles=hl.array([revel_ht.ref, revel_ht.alt]))
revel_ht = revel_ht.key_by('locus', 'alleles')
mt = mt.annotate_rows(info=mt.info.annotate(REVEL=revel_ht[mt.row_key].REVEL))

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

# annotate LOEUF from gnomAD
loeuf_v2_ht = hl.read_table(loeuf_v2_uri).key_by('transcript')
loeuf_v4_ht = hl.read_table(loeuf_v4_uri).key_by('transcript')
mt_by_transcript = mt.explode_rows(mt.vep.transcript_consequences)
mt_by_transcript = mt_by_transcript.key_rows_by(mt_by_transcript.vep.transcript_consequences.Feature)
mt_by_transcript = mt_by_transcript.annotate_rows(vep=mt_by_transcript.vep.annotate(
    transcript_consequences=mt_by_transcript.vep.transcript_consequences.annotate(
        LOEUF_v2=hl.if_else(hl.is_defined(loeuf_v2_ht[mt_by_transcript.row_key]), loeuf_v2_ht[mt_by_transcript.row_key]['oe_lof_upper'], ''),
        LOEUF_v4=hl.if_else(hl.is_defined(loeuf_v4_ht[mt_by_transcript.row_key]), loeuf_v4_ht[mt_by_transcript.row_key]['lof.oe_ci.upper'], ''))))
csq_fields_str = 'Format: ' + header['info']['CSQ']['Description'].split('Format: ')[1] + '|'.join(['', 'LOEUF_v2', 'LOEUF_v4'])

mt_by_transcript.checkpoint('loeuf_annot.mt')

mt_by_locus_and_gene = mt_by_transcript.key_rows_by('locus', 'alleles', mt_by_transcript.vep.transcript_consequences.SYMBOL)

# OPTIONAL: annotate SpliceAI scores
if (spliceAI_snv_uri!='NA') and (spliceAI_indel_uri!='NA'):
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
    csq_fields_str = csq_fields_str + '|'.join([''] + fields)
    mt_by_locus_and_gene.checkpoint('annot_spliceAI.mt')

# annotate OMIM
omim = hl.import_table(omim_uri).key_by('approvedGeneSymbol')
mt_by_gene = mt_by_locus_and_gene.key_rows_by(mt_by_locus_and_gene.vep.transcript_consequences.SYMBOL)
mt_by_gene = mt_by_gene.annotate_rows(vep=mt_by_gene.vep.annotate(
    transcript_consequences=mt_by_gene.vep.transcript_consequences.annotate(
        OMIM_MIM_number=hl.if_else(hl.is_defined(omim[mt_by_gene.row_key]), omim[mt_by_gene.row_key].mimNumber, ''),
        OMIM_inheritance_code=hl.if_else(hl.is_defined(omim[mt_by_gene.row_key]), omim[mt_by_gene.row_key].inheritance_code, ''))))
csq_fields_str = csq_fields_str + '|'.join([''] + ['OMIM_MIM_number', 'OMIM_inheritance_code'])

# OPTIONAL: annotate with gene list, if provided
if gene_list!='NA':
    genes = pd.read_csv(gene_list, sep='\t', header=None)[0].tolist()
    gene_list_name = os.path.basename(gene_list)
    mt_by_gene = mt_by_gene.annotate_rows(vep=mt_by_gene.vep.annotate(
    transcript_consequences=mt_by_gene.vep.transcript_consequences.annotate(
        gene_list=hl.if_else(hl.array(genes).contains(mt_by_gene.row_key.SYMBOL), gene_list_name, ''))))
    csq_fields_str = csq_fields_str + '|gene_list'

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