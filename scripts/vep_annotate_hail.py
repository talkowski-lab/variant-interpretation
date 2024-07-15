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
parser.add_argument('--reannotate-ac-af', dest='reannotate_ac_af', help='Whether or not AC/AF should be recalculated by Hail')
parser.add_argument('--build', dest='build', help='Genome build')
parser.add_argument('--mpc', dest='mpc_ht_uri', help='MPC scores HT')
parser.add_argument('--clinvar', dest='clinvar_vcf_uri', help='ClinVar VCF')
parser.add_argument('--omim', dest='omim_uri', help='OMIM file')
parser.add_argument('--revel', dest='revel_file', help='REVEL file')
parser.add_argument('--loeuf-v2', dest='loeuf_v2_uri', help='LOEUF scores from gnomAD v2.1.1')
parser.add_argument('--loeuf-v4', dest='loeuf_v4_uri', help='LOEUF scores from gnomAD v4.1')
parser.add_argument('--genes', dest='gene_list', help='Gene list txt file')
parser.add_argument('--project-id', dest='project_id', help='Google Project ID')

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
gcp_project = args.project_id

hl.init(min_block_size=128, 
        spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g",
        #             'spark.hadoop.fs.gs.requester.pays.mode': 'CUSTOM',
        #             'spark.hadoop.fs.gs.requester.pays.buckets': 'hail-datasets-us-central1',
        #             'spark.hadoop.fs.gs.requester.pays.project.id': gcp_project,
                    }, 
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

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
header['info']['CSQ'] = {'Description': hl.eval(mt.vep_csq_header), 'Number': '.', 'Type': 'String'}

hl.export_vcf(dataset=mt, output=vep_annotated_vcf_name, metadata=header)