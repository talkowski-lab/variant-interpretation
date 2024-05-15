#!/usr/bin/env python
# coding: utf-8

import hail as hl
import pandas as pd
import numpy as np
import sys
import os
import ast

build = 'GRCh38'
vcf_file = sys.argv[1]
ped_uri = sys.argv[2]
af_threshold = float(sys.argv[3])
cores = sys.argv[4]
mem = int(np.floor(float(sys.argv[5])))
file_ext = sys.argv[6]

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

prefix = os.path.basename(vcf_file).split(file_ext)[0]

tmp_ped = pd.read_csv(ped_uri, sep='\t')
# check ped number of columns
if len(tmp_ped) > 6:
    tmp_ped = tmp_ped.iloc[:,:6]
tmp_ped.to_csv(f"{prefix}.ped", sep='\t', index=False)

ped_uri = f"{prefix}.ped"
pedigree = hl.Pedigree.read(ped_uri, delimiter='\t')

mt = hl.import_vcf(vcf_file, force_bgz=True, array_elements_required=False, call_fields=[], reference_genome='GRCh38')

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

mt = split_multi_ssc(mt)

mt = mt.annotate_rows(info=mt.info.annotate(cohort_AC=mt.info.AC[mt.a_index - 1],
                                           cohort_AF=mt.info.AF[mt.a_index - 1]))

mt = mt.filter_entries(mt.info.cohort_AF <= af_threshold)

tm = hl.trio_matrix(mt, pedigree)

tm_denovo = tm.filter_entries((tm.proband_entry.GT.is_het()) & 
                  (tm.father_entry.GT.is_hom_ref()) & 
                  (tm.mother_entry.GT.is_hom_ref()))

tm_denovo_df = tm_denovo.entries().to_pandas()

child_format_fields = ['GT','AD','DP','DPC','GQ','PGT','PID','PL','VAF','AB']  
parent_format_fields = ['GT','AD','DP','DPC','GQ','VAF','AB']  
rename_cols = {f"mother_entry.{field}": f"{field}_mother" for field in parent_format_fields if f"mother_entry.{field}" in tm_denovo_df.columns} |\
    {f"father_entry.{field}": f"{field}_father" for field in parent_format_fields if f"father_entry.{field}" in tm_denovo_df.columns} |\
    {f"proband_entry.{field}": f"{field}_sample" for field in child_format_fields if f"proband_entry.{field}" in tm_denovo_df.columns} |\
    {'qual': 'QUAL', 'proband.s': 'SAMPLE', 'filters': 'FILTER'}
tm_denovo_df = tm_denovo_df.rename(rename_cols, axis=1)

tm_denovo_df['CHROM'] = tm_denovo_df.locus.astype(str).str.split(':').str[0]
tm_denovo_df['POS'] = tm_denovo_df.locus.astype(str).str.split(':').str[1].astype(int)
tm_denovo_df['REF'] = tm_denovo_df.alleles.str[0]
tm_denovo_df['ALT'] = tm_denovo_df.alleles.str[1]
tm_denovo_df['LEN'] = abs(tm_denovo_df.REF.str.len()-tm_denovo_df.ALT.str.len())
tm_denovo_df['TYPE'] =np.where(tm_denovo_df.LEN==0, 'SNV', 'Indel')

tm_denovo_df.columns = tm_denovo_df.columns.str.replace('info.', '')

for info_cat in ['AC', 'AF', 'MLEAC', 'MLEAF']:
    if info_cat in tm_denovo_df.columns:
            tm_denovo_df[info_cat] = tm_denovo_df[info_cat].str[0]

# 'POLYX' -- added after downsampling
info_cols = ['END','AC','AF','AN','BaseQRankSum','ClippingRankSum','DP','FS','MLEAC','MLEAF','MQ','MQRankSum','QD','ReadPosRankSum','SOR','VQSLOD','cohort_AC', 'cohort_AF', 'CSQ']
info_cols = list(np.intersect1d(info_cols, list(mt.info.keys())))
cols_to_keep = ['CHROM', 'POS', 'REF', 'ALT', 'LEN', 'TYPE', 'ID', 'VarKey'] + info_cols + list(rename_cols.values())

tm_denovo_df[cols_to_keep].to_csv(f"{os.path.basename(vcf_file).split(file_ext)[0]}_denovo_GT_AF_filter.tsv.gz", sep='\t', index=False)