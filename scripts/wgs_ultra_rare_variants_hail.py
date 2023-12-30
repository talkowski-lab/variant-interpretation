#!/usr/bin/env python
# coding: utf-8

# adapted from wgs_ultra_rare_variants_hail.ipynb in Terra Analysis

import hail as hl
import pandas as pd
import numpy as np
import sys

build = 'GRCh38'
lcr_uri = sys.argv[1]
ped_uri = sys.argv[2]
meta_uri = sys.argv[3]
trio_uri = sys.argv[4]
vcf_file = sys.argv[5]
cohort_prefix = sys.argv[6]
cores = sys.argv[7]
mem = int(np.floor(float(sys.argv[8])))
ac_threshold = int(sys.argv[9])
af_threshold = float(sys.argv[10])

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.memory": f"{mem}g"})

pedigree = hl.Pedigree.read(ped_uri)

trio_df = pd.read_csv(trio_uri, sep='\t')

mt = hl.import_vcf(vcf_file, force_bgz=True, reference_genome='GRCh38')

#split-multi
def split_multi_ssc(mt):
    mt = mt.annotate_rows(num_alleles = mt.alleles.size() ) # Add number of alleles at site before split
    # Now split
    sm = hl.split_multi(mt)
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
# annotate cohort ac to INFO field (after splitting multiallelic)
mt = mt.annotate_rows(info=mt.info.annotate(cohort_AC=mt.info.AC[mt.a_index - 1],
                                           cohort_AF=mt.info.AF[mt.a_index - 1]))
mt = mt.annotate_rows(ID=hl.variant_str(mt.row_key))

# try filtering before
mt_filtered = mt
mt_filtered = mt_filtered.filter_rows(mt_filtered.info.cohort_AC < 20 , keep=True)
# filter low complexity regions
try:
    lcr = hl.import_bed(lcr_uri, reference_genome=build)
except Exception as e:
    lcr = hl.import_bed(lcr_uri, reference_genome=build, force_bgz=True)
mt_filtered = mt_filtered.filter_rows(hl.is_defined(lcr[mt_filtered.locus]), keep=False)
# select samples after qc
meta = hl.import_table(meta_uri, types={'SampleID': hl.tstr, 'FamID': hl.tstr, 'Role': hl.tstr})
meta = meta.annotate(s=meta.SampleID).key_by('s')
mt_filtered = mt_filtered.annotate_cols(pheno=meta[mt_filtered.s])
mt_filtered = mt_filtered.filter_cols(mt_filtered.pheno.Role != '', keep = True)

# # filter for VQSR - PASS for SNVs only
mt_filtered = mt_filtered.filter_rows((hl.is_snp(mt_filtered.alleles[0], mt_filtered.alleles[1]) &(mt_filtered.filters.size() == 0))
                   | hl.is_indel(mt_filtered.alleles[0], mt_filtered.alleles[1]))
# filter on depth
mt_filtered = mt_filtered.filter_entries( (mt_filtered.DP < 10) | (mt_filtered.DP > 200), keep = False) 
# row (variant INFO) level filters - GATK recommendations for short variants
dn_snv_cond_row = (hl.is_snp(mt_filtered.alleles[0], mt_filtered.alleles[1])
                    & (mt_filtered.qual >= 150)
                    & (mt_filtered.info.SOR <= 2.5)
                    & (mt_filtered.info.ReadPosRankSum >= -1.4)
                    & (mt_filtered.info.QD >= 3.0)
                    & (mt_filtered.info.MQ >= 50))
dn_indel_cond_row = (hl.is_indel(mt_filtered.alleles[0], mt_filtered.alleles[1])
                      & (mt_filtered.qual >= 150)
                      & (mt_filtered.info.SOR <= 3)
                      & (mt_filtered.info.ReadPosRankSum >= -1.7)
                      & (mt_filtered.info.QD >= 4.0)
                      & (mt_filtered.info.MQ >= 50))
mt_filtered = mt_filtered.filter_rows(dn_snv_cond_row | dn_indel_cond_row, keep = True)

# GQ mean filters
mt_filtered = hl.variant_qc(mt_filtered)
mt_filtered = mt_filtered.filter_rows(mt_filtered.variant_qc.gq_stats.mean >= 50, keep = True)
# clean-up: remove AC = 0 loci
mt_filtered = hl.variant_qc(mt_filtered)
mt_filtered = mt_filtered.filter_rows(mt_filtered.variant_qc.AC[1] > 0, keep = True)
# write to output vcf
mt_filtered = mt_filtered.drop('variant_qc')

# useful table here: https://hail.is/docs/0.2/methods/genetics.html#hail.methods.transmission_disequilibrium_test
# t=1 u=0
tdt_table_filtered = hl.transmission_disequilibrium_test(mt_filtered, pedigree)

mt_filtered_rare = mt_filtered.filter_rows((mt_filtered.info.cohort_AC<=ac_threshold)|(mt_filtered.info.cohort_AF<=af_threshold))
tdt_table_filtered_rare = tdt_table_filtered.semi_join(mt_filtered_rare.rows())

ultra_rare_vars_table = tdt_table_filtered_rare.filter((tdt_table_filtered_rare.t==1) & (tdt_table_filtered_rare.u==0))

trio_mat = hl.trio_matrix(mt_filtered, pedigree, complete_trios=True).semi_join_rows(ultra_rare_vars_table)

ultra_rare_vars_df = trio_mat.filter_entries(trio_mat.proband_entry.GT.is_het()).entries().to_pandas()

ultra_rare_vars_df = ultra_rare_vars_df[ultra_rare_vars_df['proband.s'].isin(trio_df.SampleID)]

ultra_rare_vars_df['father_entry.VAF'] = ultra_rare_vars_df['father_entry.AD'].str[1] / (ultra_rare_vars_df['father_entry.AD'].str[0]+ultra_rare_vars_df['father_entry.AD'].str[1])
ultra_rare_vars_df['mother_entry.VAF'] = ultra_rare_vars_df['mother_entry.AD'].str[1] / (ultra_rare_vars_df['mother_entry.AD'].str[0]+ultra_rare_vars_df['mother_entry.AD'].str[1])
ultra_rare_vars_df['proband_entry.VAF'] = ultra_rare_vars_df['proband_entry.AD'].str[1] / (ultra_rare_vars_df['proband_entry.AD'].str[0]+ultra_rare_vars_df['proband_entry.AD'].str[1])

child_format_fields = ['GT','AD','DP','GQ','PGT','PID','PL','VAF']  #'AB'
parent_format_fields = ['GT','AD','DP','GQ','VAF']  #'AB'
rename_cols = {f"mother_entry.{field}": f"{field}_mother" for field in parent_format_fields} |\
    {f"father_entry.{field}": f"{field}_father" for field in parent_format_fields} |\
    {f"proband_entry.{field}": f"{field}_sample" for field in child_format_fields} |\
    {'qual': 'QUAL', 'proband.s': 'SAMPLE'}
ultra_rare_vars_df = ultra_rare_vars_df.rename(rename_cols, axis=1)

ultra_rare_vars_df['CHROM'] = ultra_rare_vars_df.locus.astype(str).str.split(':').str[0]
ultra_rare_vars_df['POS'] = ultra_rare_vars_df.locus.astype(str).str.split(':').str[1].astype(int)
ultra_rare_vars_df['REF'] = ultra_rare_vars_df.alleles.str[0]
ultra_rare_vars_df['ALT'] = ultra_rare_vars_df.alleles.str[1]
ultra_rare_vars_df['LEN'] = abs(ultra_rare_vars_df.REF.str.len()-ultra_rare_vars_df.ALT.str.len())
ultra_rare_vars_df['TYPE'] =np.where(ultra_rare_vars_df.LEN==0, 'SNV', 'Indel')

ultra_rare_vars_df.columns = ultra_rare_vars_df.columns.str.replace('info.', '')

# 'POLYX'
info_cols = ['END','AC','AF','AN','BaseQRankSum','ClippingRankSum','DP','FS','MLEAC','MLEAF','MQ','MQRankSum','QD','ReadPosRankSum','SOR','VQSLOD','cohort_AC', 'cohort_AF', 'CSQ']
cols_to_keep = ['CHROM', 'POS', 'REF', 'ALT', 'LEN', 'TYPE', 'ID'] + info_cols + list(rename_cols.values())

ultra_rare_vars_df[cols_to_keep].to_csv(f"{cohort_prefix}_ultra_rare_variants.tsv", sep='\t', index=False)