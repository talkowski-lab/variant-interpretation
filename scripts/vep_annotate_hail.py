from pyspark.sql import SparkSession
import hail as hl
import numpy as np
import pandas as pd
import sys
import ast
import os
import argparse

parser = argparse.ArgumentParser(description='Parse arguments')
parser.add_argument('-i', dest='vcf_file', help='Input VCF file')
parser.add_argument('-o', dest='vep_annotated_vcf_name', help='Output filename')
parser.add_argument('--cores', dest='cores', help='CPU cores')
parser.add_argument('--mem', dest='mem', help='Memory')
parser.add_argument('--reannotate-ac-af', dest='reannotate_ac_af', help='Whether or not AC/AF should be recalculated by Hail')
parser.add_argument('--build', dest='build', help='Genome build')
parser.add_argument('--mpc', dest='mpc_ht_uri', help='MPC scores HT')
parser.add_argument('--clinvar', dest='clinvar_vcf_uri', help='ClinVar VCF')
parser.add_argument('--omim', dest='omim_uri', help='OMIM file')
parser.add_argument('--revel', dest='revel_file', help='REVEL file')
parser.add_argument('--loeuf-v2', dest='loeuf_v2_uri', help='LOEUF scores from gnomAD v2.1.1')
parser.add_argument('--loeuf-v4', dest='loeuf_v4_uri', help='LOEUF scores from gnomAD v4.1')
parser.add_argument('--genes', dest='gene_list', help='Gene list txt file')

args = parser.parse_args()

vcf_file = args.vcf_file
vep_annotated_vcf_name = args.vep_annotated_vcf_name
cores = args.cores  # string
mem = int(np.floor(float(args.mem)))
reannotate_ac_af = ast.literal_eval(args.reannotate_ac_af.capitalize())
build = args.build
mpc_ht_uri = args.mpc_ht_uri
clinvar_vcf_uri = args.clinvar_vcf_uri
omim_uri = args.omim_uri
revel_file = args.revel_file
loeuf_v2_uri = args.loeuf_v2_uri
loeuf_v4_uri = args.loeuf_v4_uri
gene_list = args.gene_list

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

#split-multi
def split_multi_ssc(mt):
    mt = mt.annotate_rows(num_alleles = mt.alleles.size() ) # Add number of alleles at site before split
    # only split variants that aren't already split
    bi = mt.filter_rows(hl.len(mt.alleles) == 2)
    bi = bi.annotate_rows(a_index=1, was_split=False, old_locus=bi.locus, old_alleles=bi.alleles)
    multi = mt.filter_rows(hl.len(mt.alleles) > 2)
    # Now split
    split = hl.split_multi(multi, permit_shuffle=True)
    sm = split.union_rows(bi)
    # sm = hl.split_multi(mt, permit_shuffle=True)
    pl = hl.or_missing(hl.is_defined(sm.PL),
                      (hl.range(0, 3).map(lambda i: hl.min(hl.range(0, hl.len(sm.PL))
       .filter(lambda j: hl.downcode(hl.unphased_diploid_gt_index_call(j), sm.a_index) == hl.unphased_diploid_gt_index_call(i))
       .map(lambda j: sm.PL[j])))))
    split_ds = sm.annotate_entries(GT = hl.downcode(sm.GT, sm.a_index),
                                   AD = hl.or_missing(hl.is_defined(sm.AD), [hl.sum(sm.AD) - sm.AD[sm.a_index], sm.AD[sm.a_index]]),
                                   PL = pl) 
        #GQ = hl.cond(hl.is_defined(pl[0]) & hl.is_defined(pl[1]) & hl.is_defined(pl[2]), hl.gq_from_pl(pl), sm.GQ) )
    mt = split_ds.drop('old_locus', 'old_alleles')
    return mt

header = hl.get_vcf_metadata(vcf_file) 
mt = hl.import_vcf(vcf_file, force_bgz=True, array_elements_required=False, call_fields=[], reference_genome=build)

# mt = mt.distinct_by_row()
if 'num_alleles' not in list(mt.row_value.keys()):
    mt = split_multi_ssc(mt)
    # mt = mt.distinct_by_row()

try:
    # for haploid (e.g. chrY)
    mt = mt.annotate_entries(
        GT = hl.if_else(
                mt.GT.ploidy == 1, 
                hl.call(mt.GT[0], mt.GT[0]),
                mt.GT)
    )
except:
    pass

# annotate cohort AC to INFO field (after splitting multiallelic)
mt = mt.annotate_rows(info=mt.info.annotate(cohort_AC=mt.info.AC[mt.a_index - 1],
                                            cohort_AF=mt.info.AF[mt.a_index - 1]))

# reannotate
if (reannotate_ac_af):
    mt = hl.variant_qc(mt)
    mt = mt.annotate_rows(info=mt.info.annotate(cohort_AC=mt.variant_qc.AC[1],
                                      cohort_AF=mt.variant_qc.AF[1]))
    mt = mt.drop('variant_qc')

# for VCFs with AS_VQSLOD and missing VQSLOD
all_as_fields = [col for col in list(mt.info) if 'AS_' in col]
for field in all_as_fields:
    normal_field = field.split('_')[1]
    n_missing_as = mt.filter_rows(hl.is_missing(getattr(mt.info, field))).count_rows()
    if normal_field not in list(mt.info):
        continue
    n_missing = mt.filter_rows(hl.is_missing(getattr(mt.info, normal_field))).count_rows()
    if (n_missing_as < n_missing):
        mt = mt.annotate_rows(info=mt.info.annotate(**{normal_field: getattr(mt.info, field)[mt.a_index - 1]}))    

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

# run VEP
mt = hl.vep(mt, config='vep_config.json', csq=True, tolerate_parse_error=True)
mt = mt.annotate_rows(info = mt.info.annotate(CSQ=mt.vep))

# annotate REVEL
revel_ht = hl.import_table(revel_file, force_bgz=True)
revel_ht = revel_ht.annotate(chr='chr'+revel_ht['#chr']) 
build_chr = 'chr' if build=='GRCh38' else '#chr'
build_pos = 'grch38_pos' if build=='GRCh38' else 'hg19_pos'
revel_ht = revel_ht.annotate(locus=hl.locus(revel_ht[build_chr], hl.int(revel_ht[build_pos]), build),
                 alleles=hl.array([revel_ht.ref, revel_ht.alt]))
revel_ht = revel_ht.key_by('locus', 'alleles')
mt = mt.annotate_rows(info=mt.info.annotate(REVEL=revel_ht[mt.row_key].REVEL))

# annotate OMIM
csq_columns = hl.eval(mt.vep_csq_header).split('Format: ')[1].split('|')
mt = mt.annotate_rows(vep=mt.info)
transcript_consequences = mt.vep.CSQ.map(lambda x: x.split('\|'))

transcript_consequences_strs = transcript_consequences.map(lambda x: hl.if_else(hl.len(x)>1, hl.struct(**
                                                       {col: x[i] if col!='Consequence' else x[i].split('&')  
                                                        for i, col in enumerate(csq_columns)}), 
                                                        hl.struct(**{col: hl.missing('str') if col!='Consequence' else hl.array([hl.missing('str')])  
                                                        for i, col in enumerate(csq_columns)})))

mt = mt.annotate_rows(vep=mt.vep.annotate(transcript_consequences=transcript_consequences_strs))
mt = mt.annotate_rows(vep=mt.vep.select('transcript_consequences'))

omim = hl.import_table(omim_uri).key_by('approvedGeneSymbol')

mt_by_gene = mt.explode_rows(mt.vep.transcript_consequences)
mt_by_gene = mt_by_gene.key_rows_by(mt_by_gene.vep.transcript_consequences.SYMBOL)
mt_by_gene = mt_by_gene.annotate_rows(vep=mt_by_gene.vep.annotate(
    transcript_consequences=mt_by_gene.vep.transcript_consequences.annotate(
        OMIM_MIM_number=hl.if_else(hl.is_defined(omim[mt_by_gene.row_key]), omim[mt_by_gene.row_key].mimNumber, ''),
        OMIM_inheritance_code=hl.if_else(hl.is_defined(omim[mt_by_gene.row_key]), omim[mt_by_gene.row_key].inheritance_code, ''))))

# annotate LOEUF from gnomAD
loeuf_v4_ht = hl.import_table(loeuf_v4_uri, force_bgz=loeuf_v4_uri.split('.')[-1] in ['bgz','gz']).key_by('transcript')
loeuf_v2_ht = hl.import_table(loeuf_v2_uri, force_bgz=loeuf_v2_uri.split('.')[-1] in ['bgz','gz']).key_by('transcript')
mt_by_transcript = mt_by_gene.key_rows_by(mt_by_gene.vep.transcript_consequences.Feature)
mt_by_transcript = mt_by_transcript.annotate_rows(vep=mt_by_transcript.vep.annotate(
    transcript_consequences=mt_by_transcript.vep.transcript_consequences.annotate(
        LOEUF_v2=hl.if_else(hl.is_defined(loeuf_v2_ht[mt_by_transcript.row_key]), loeuf_v2_ht[mt_by_transcript.row_key]['oe_lof_upper'], ''),
        LOEUF_v4=hl.if_else(hl.is_defined(loeuf_v4_ht[mt_by_transcript.row_key]), loeuf_v4_ht[mt_by_transcript.row_key]['lof.oe_ci.upper'], ''))))

csq_fields_str = hl.eval(mt.vep_csq_header) + '|'.join(['', 'OMIM_MIM_number', 'OMIM_inheritance_code', 'LOEUF_v2', 'LOEUF_v4'])

# annotate with gene list, if provided
if gene_list.split('.')[-1] == 'txt':
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

hl.export_vcf(dataset=mt, output=vep_annotated_vcf_name, metadata=header)