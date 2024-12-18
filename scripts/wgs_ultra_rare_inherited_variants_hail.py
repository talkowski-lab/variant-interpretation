#!/usr/bin/env python
# coding: utf-8

# adapted from wgs_ultra_rare_variants_hail.ipynb in Terra Analysis

import hail as hl
import pandas as pd
import numpy as np
import sys
import ast
import warnings
import os

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
csq_af_threshold = float(sys.argv[11])
gq_het_threshold = float(sys.argv[12])
gq_hom_ref_threshold = float(sys.argv[13])
qual_threshold = int(sys.argv[14])
sor_threshold_indel = float(sys.argv[15])
sor_threshold_snv = float(sys.argv[16])
readposranksum_threshold_indel = float(sys.argv[17])
readposranksum_threshold_snv = float(sys.argv[18])
qd_threshold_indel = float(sys.argv[19])
qd_threshold_snv = float(sys.argv[20])
mq_threshold = float(sys.argv[21])
build = sys.argv[22]

hl.init(min_block_size=128, 
        local=f"local[*]", 
        spark_conf={
                    "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.speculation": 'true'
                    }, 
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

trio_df = pd.read_csv(trio_uri, dtype=str, sep='\t')

prefix = os.path.basename(vcf_file).split('.vcf')[0]
mt = hl.import_vcf(vcf_file, force_bgz=True, array_elements_required=False, call_fields=[], reference_genome=build)

tmp_ped = pd.read_csv(ped_uri, sep='\t')
# check ped number of columns
if len(tmp_ped) > 6:
    tmp_ped = tmp_ped.iloc[:,:6]
    
# subset ped to samples in mt
samps = mt.s.collect()
tmp_ped = tmp_ped[tmp_ped.iloc[:,1].isin(samps)]  # sample_id
tmp_ped.columns = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
tmp_ped = tmp_ped.drop_duplicates('sample_id')    
tmp_ped = tmp_ped.replace({np.nan: 0})

tmp_ped.to_csv(f"{prefix}.ped", sep='\t', index=False)
pedigree = hl.Pedigree.read(f"{prefix}.ped", delimiter='\t')

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
    if 'PL' in list(mt.entry.keys()):
        pl = hl.or_missing(hl.is_defined(sm.PL),
                        (hl.range(0, 3).map(lambda i: hl.min(hl.range(0, hl.len(sm.PL))
        .filter(lambda j: hl.downcode(hl.unphased_diploid_gt_index_call(j), sm.a_index) == hl.unphased_diploid_gt_index_call(i))
        .map(lambda j: sm.PL[j])))))
        sm = sm.annotate_entries(PL = pl)
    split_ds = sm.annotate_entries(GT = hl.downcode(sm.GT, sm.a_index),
                                   AD = hl.or_missing(hl.is_defined(sm.AD), [hl.sum(sm.AD) - sm.AD[sm.a_index], sm.AD[sm.a_index]])
                                   ) 
        #GQ = hl.cond(hl.is_defined(pl[0]) & hl.is_defined(pl[1]) & hl.is_defined(pl[2]), hl.gq_from_pl(pl), sm.GQ) )
    mt = split_ds.drop('old_locus', 'old_alleles')
    return mt

mt = split_multi_ssc(mt)
# annotate cohort ac to INFO field (after splitting multiallelic)
if 'cohort_AC' not in list(mt.info):
    mt = mt.annotate_rows(info=mt.info.annotate(cohort_AC=mt.info.AC[mt.a_index - 1]))
if 'cohort_AF' not in list(mt.info):
    mt = mt.annotate_rows(info=mt.info.annotate(cohort_AF=mt.info.AF[mt.a_index - 1]))
mt = mt.annotate_rows(ID=hl.variant_str(mt.row_key))

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

# # filter for VQSR - PASS 
mt_filtered = mt_filtered.filter_rows((mt_filtered.filters.size() == 0))
# filter on depth
mt_filtered = mt_filtered.annotate_entries(DPC=hl.sum(mt_filtered.AD),
                                           AB=mt_filtered.AD[1]/hl.sum(mt_filtered.AD),
                                           VAF=mt_filtered.AD[1]/mt_filtered.DP)
mt_filtered = mt_filtered.filter_entries( (mt_filtered.DPC < 10) | (mt_filtered.DPC > 200), keep = False) 

# row (variant INFO) level filters - GATK recommendations for short variants
dn_snv_cond_row = (hl.is_snp(mt_filtered.alleles[0], mt_filtered.alleles[1])
                    & (mt_filtered.qual >= qual_threshold)
                    & (mt_filtered.info.SOR <= sor_threshold_snv)
                    & (mt_filtered.info.ReadPosRankSum >= readposranksum_threshold_snv)
                    & (mt_filtered.info.QD >= qd_threshold_snv)
                    & (mt_filtered.info.MQ >= mq_threshold))
dn_indel_cond_row = (hl.is_indel(mt_filtered.alleles[0], mt_filtered.alleles[1])
                    & (mt_filtered.qual >= qual_threshold)
                    & (mt_filtered.info.SOR <= sor_threshold_indel)
                    & (mt_filtered.info.ReadPosRankSum >= readposranksum_threshold_indel)
                    & (mt_filtered.info.QD >= qd_threshold_indel)
                    & (mt_filtered.info.MQ >= mq_threshold))
mt_filtered = mt_filtered.filter_rows(dn_snv_cond_row | dn_indel_cond_row, keep = True)

# clean-up: remove AC = 0 loci
mt_filtered = hl.variant_qc(mt_filtered)
mt_filtered = mt_filtered.filter_rows(mt_filtered.variant_qc.AC[1] > 0, keep = True)
mt_filtered = mt_filtered.drop('variant_qc')

# useful table here: https://hail.is/docs/0.2/methods/genetics.html#hail.methods.transmission_disequilibrium_test
# t=1 u=0
tdt_table_filtered = hl.transmission_disequilibrium_test(mt_filtered, pedigree)

mt_filtered_rare = mt_filtered.filter_rows((mt_filtered.info.cohort_AC<=ac_threshold)|(mt_filtered.info.cohort_AF<=af_threshold))
tdt_table_filtered_rare = tdt_table_filtered.semi_join(mt_filtered_rare.rows())

ultra_rare_vars_table = tdt_table_filtered_rare.filter((tdt_table_filtered_rare.t==1) & (tdt_table_filtered_rare.u==0))

trio_mat = hl.trio_matrix(mt_filtered_rare, pedigree, complete_trios=True).semi_join_rows(ultra_rare_vars_table)


trio_mat = trio_mat.filter_entries((trio_mat.proband_entry.GT.is_het() 
                                            & (trio_mat.father_entry.GT.is_het() 
                                               | trio_mat.mother_entry.GT.is_het())))

# GQ sample/parent filters
trio_mat = trio_mat.filter_entries(hl.if_else(trio_mat.father_entry.GT.is_het(), trio_mat.father_entry.GQ>=gq_het_threshold,
                                              trio_mat.father_entry.GQ>=gq_hom_ref_threshold)
                                    & hl.if_else(trio_mat.mother_entry.GT.is_het(), trio_mat.mother_entry.GQ>=gq_het_threshold,
                                                trio_mat.mother_entry.GQ>=gq_hom_ref_threshold) 
                                    & (trio_mat.proband_entry.GQ>=gq_het_threshold))

# AB sample/parent filters
trio_mat = trio_mat.filter_entries(hl.if_else(trio_mat.father_entry.GT.is_het(), (trio_mat.father_entry.AB>=0.2) & (trio_mat.father_entry.AB<=0.8),
                                              trio_mat.father_entry.AB<=0.05)
                                    & hl.if_else(trio_mat.mother_entry.GT.is_het(), (trio_mat.mother_entry.AB>=0.2) & (trio_mat.mother_entry.AB<=0.8),
                                                trio_mat.mother_entry.AB<=0.05) 
                                    & ((trio_mat.proband_entry.AB>=0.2) & (trio_mat.proband_entry.AB<=0.8)))

# DP filter for indels
trio_mat = trio_mat.filter_entries(((hl.is_indel(trio_mat.alleles[0], trio_mat.alleles[1])) 
                                        & (trio_mat.mother_entry.DPC>=16) 
                                        & (trio_mat.father_entry.DPC>=16)) 
                                   | (hl.is_snp(trio_mat.alleles[0], trio_mat.alleles[1])))

# get random sample's GQ
mt_ultra_rare = mt.semi_join_rows(trio_mat.rows())
mt_ultra_rare = mt_ultra_rare.filter_entries((mt_ultra_rare.GT.is_hom_ref()) & (mt_ultra_rare.GQ>=gq_hom_ref_threshold))

# clean-up: remove AC = 0 loci
mt_ultra_rare = hl.variant_qc(mt_ultra_rare)
mt_ultra_rare = mt_ultra_rare.filter_rows(mt_ultra_rare.variant_qc.AC[0] > 0, keep = True) 
mt_ultra_rare = mt_ultra_rare.drop('variant_qc')

mt_ultra_rare = mt_ultra_rare.annotate_rows(all_GQs=hl.array(hl.agg.collect_as_set(mt_ultra_rare.GQ)))
mt_ultra_rare = mt_ultra_rare.annotate_rows(GQ_random=hl.shuffle(mt_ultra_rare.all_GQs)[0])

trio_mat = trio_mat.annotate_rows(GQ_random=mt_ultra_rare.rows()[trio_mat.row_key].GQ_random)
trio_mat = trio_mat.filter_rows(hl.is_defined(trio_mat.GQ_random))

ultra_rare_vars_df = trio_mat.entries().to_pandas()

ultra_rare_vars_df = ultra_rare_vars_df[ultra_rare_vars_df['proband.s'].isin(trio_df.SampleID)]

# ultra_rare_vars_df['father_entry.VAF'] = ultra_rare_vars_df['father_entry.AD'].str[1] / ultra_rare_vars_df['father_entry.DP']
# ultra_rare_vars_df['mother_entry.VAF'] = ultra_rare_vars_df['mother_entry.AD'].str[1] / ultra_rare_vars_df['mother_entry.DP']
# ultra_rare_vars_df['proband_entry.VAF'] = ultra_rare_vars_df['proband_entry.AD'].str[1] / ultra_rare_vars_df['proband_entry.DP']

child_format_fields = ['GT','AD','DP','DPC','GQ','PGT','PID','PL','VAF','AB']  
parent_format_fields = ['GT','AD','DP','DPC','GQ','VAF','AB']  
rename_cols = {f"mother_entry.{field}": f"{field}_mother" for field in parent_format_fields if f"mother_entry.{field}" in ultra_rare_vars_df.columns} |\
    {f"father_entry.{field}": f"{field}_father" for field in parent_format_fields if f"father_entry.{field}" in ultra_rare_vars_df.columns} |\
    {f"proband_entry.{field}": f"{field}_sample" for field in child_format_fields if f"proband_entry.{field}" in ultra_rare_vars_df.columns} |\
    {'qual': 'QUAL', 'proband.s': 'SAMPLE', 'filters': 'FILTER'}
ultra_rare_vars_df = ultra_rare_vars_df.rename(rename_cols, axis=1)

ultra_rare_vars_df['CHROM'] = ultra_rare_vars_df.locus.astype(str).str.split(':').str[0]
ultra_rare_vars_df['POS'] = ultra_rare_vars_df.locus.astype(str).str.split(':').str[1].astype(int)
ultra_rare_vars_df['REF'] = ultra_rare_vars_df.alleles.str[0]
ultra_rare_vars_df['ALT'] = ultra_rare_vars_df.alleles.str[1]
ultra_rare_vars_df['LEN'] = abs(ultra_rare_vars_df.REF.str.len()-ultra_rare_vars_df.ALT.str.len())
ultra_rare_vars_df['TYPE'] =np.where(ultra_rare_vars_df.LEN==0, 'SNV', 'Indel')

ultra_rare_vars_df.columns = ultra_rare_vars_df.columns.str.replace('info.', '')

# for info_cat in ['AC', 'AF', 'MLEAC', 'MLEAF']:
#     if info_cat in ultra_rare_vars_df.columns:
#             ultra_rare_vars_df[info_cat] = ultra_rare_vars_df[info_cat].str[0]

# 'POLYX' -- added after downsampling
info_cols = ['END','AC','AF','AN','BaseQRankSum','ClippingRankSum','DP','FS','MLEAC','MLEAF','MQ','MQRankSum','QD','ReadPosRankSum','SOR','VQSLOD','cohort_AC', 'cohort_AF', 'CSQ']
info_cols = list(np.intersect1d(info_cols, list(mt_filtered_rare.info.keys())))
cols_to_keep = ['CHROM', 'POS', 'REF', 'ALT', 'LEN', 'TYPE', 'ID', 'VarKey', 'GQ_random'] + info_cols + list(rename_cols.values())
# cols_to_keep = ['CHROM', 'POS', 'REF', 'ALT', 'LEN', 'TYPE', 'ID', 'VarKey'] + info_cols + list(rename_cols.values())

# CSQ AF threshold
def get_gnomAD_AF(csq, col_num):
    if type(csq)==float:
        return np.nan
    csqs = []
    for ind_csq in csq:
        af = ind_csq.split('|')[col_num]
        if af != '':
            csqs.append(af)
    csqs = list(set(csqs))
    if len(csqs)==0:
        return np.nan
    return csqs[0]

ultra_rare_vars_df['CSQ'] = ultra_rare_vars_df.CSQ.replace({'.':np.nan, None: np.nan})#.str.split(',')

try:
    header = hl.get_vcf_metadata(vcf_file)
    csq_columns = header['info']['CSQ']['Description'].split('Format: ')[1].split('|')

    for gnomad_af_str in ['gnomADe_AF', 'gnomADg_AF']:
        ultra_rare_vars_df[gnomad_af_str] = ultra_rare_vars_df.CSQ.apply(get_gnomAD_AF, col_num=csq_columns.index(gnomad_af_str)).astype(float)
    
    ultra_rare_vars_df['gnomAD_max_AF'] = ultra_rare_vars_df[['gnomADe_AF', 'gnomADg_AF']].max(axis=1)
    ultra_rare_vars_df = ultra_rare_vars_df[ultra_rare_vars_df['gnomAD_max_AF'].replace({np.nan: 0})<=csq_af_threshold]
    cols_to_keep = cols_to_keep + ['gnomADe_AF', 'gnomADg_AF', 'gnomAD_max_AF']

except Exception as e:
    print(str(e))
    # pass

try:
    ultra_rare_vars_df['VarKey'] = ultra_rare_vars_df[['ID', 'SAMPLE']].astype(str).agg(':'.join, axis=1)
except Exception as e:
    print(str(e))
    ultra_rare_vars_df['VarKey'] = np.nan
ultra_rare_vars_df[cols_to_keep].to_csv(f"{cohort_prefix}_ultra_rare_inherited_variants.tsv.gz", sep='\t', index=False)
