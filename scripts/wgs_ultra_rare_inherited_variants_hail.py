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
csq_af_threshold = float(sys.argv[11])
gq_het_threshold = float(sys.argv[12])
gq_hom_ref_threshold = float(sys.argv[13])

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

trio_df = pd.read_csv(trio_uri, dtype=str, sep='\t')

prefix = os.path.basename(vcf_file).split('.vcf')[0]
mt = hl.import_vcf(vcf_file, force_bgz=True, array_elements_required=False, call_fields=[], reference_genome='GRCh38')

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

# only keep mendelian errors 
# mendelian errors code == 2 get parents hom_ref and children het loci
pedigree = hl.Pedigree.read(f"{prefix}.ped", delimiter='\t')

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

# # filter for VQSR - PASS 
mt_filtered = mt_filtered.filter_rows((mt_filtered.filters.size() == 0))
# filter on depth
mt_filtered = mt_filtered.annotate_entries(DPC=hl.sum(mt_filtered.AD),
                                           AB=mt_filtered.AD[1]/hl.sum(mt_filtered.AD),
                                           VAF=mt_filtered.AD[1]/mt_filtered.DP)
mt_filtered = mt_filtered.filter_entries( (mt_filtered.DPC < 10) | (mt_filtered.DPC > 200), keep = False) 

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

for info_cat in ['AC', 'AF', 'MLEAC', 'MLEAF']:
    if info_cat in ultra_rare_vars_df.columns:
            ultra_rare_vars_df[info_cat] = ultra_rare_vars_df[info_cat].str[0]

# 'POLYX' -- added after downsampling
info_cols = ['END','AC','AF','AN','BaseQRankSum','ClippingRankSum','DP','FS','MLEAC','MLEAF','MQ','MQRankSum','QD','ReadPosRankSum','SOR','VQSLOD','cohort_AC', 'cohort_AF', 'CSQ']
info_cols = list(np.intersect1d(info_cols, list(mt_filtered_rare.info.keys())))
cols_to_keep = ['CHROM', 'POS', 'REF', 'ALT', 'LEN', 'TYPE', 'ID', 'VarKey', 'GQ_random'] + info_cols + list(rename_cols.values())

# CSQ AF threshold
csq_columns_less = ['Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 
                    'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 
                    'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'ALLELE_NUM', 
                    'DISTANCE', 'STRAND', 'FLAGS', 'VARIANT_CLASS', 'MINIMISED', 'SYMBOL_SOURCE', 
                    'HGNC_ID', 'CANONICAL', 'TSL', 'APPRIS', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 
                    'UNIPARC', 'GENE_PHENO', 'SIFT', 'PolyPhen', 'DOMAINS', 'miRNA', 'HGVS_OFFSET', 
                    'AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'AA_AF', 'EA_AF', 'gnomAD_AF', 
                    'gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 
                    'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF', 'MAX_AF', 'MAX_AF_POPS', 'CLIN_SIG', 
                    'SOMATIC', 'PHENO', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 
                    'LOEUF', 'LoF', 'LoF_filter', 'LoF_flags', 'LoF_info']

csq_columns_more = ["Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature",
                   "BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position",
                   "Protein_position","Amino_acids","Codons","Existing_variation","ALLELE_NUM",
                   "DISTANCE","STRAND","FLAGS","VARIANT_CLASS","MINIMISED","SYMBOL_SOURCE",
                   "HGNC_ID","CANONICAL","MANE_SELECT","MANE_PLUS_CLINICAL","TSL","APPRIS",
                   "CCDS","ENSP","SWISSPROT","TREMBL","UNIPARC","UNIPROT_ISOFORM","GENE_PHENO",
                   "SIFT","PolyPhen","DOMAINS","miRNA","HGVS_OFFSET","AF","AFR_AF","AMR_AF",
                   "EAS_AF","EUR_AF","SAS_AF","gnomADe_AF","gnomADe_AFR_AF","gnomADe_AMR_AF",
                   "gnomADe_ASJ_AF","gnomADe_EAS_AF","gnomADe_FIN_AF","gnomADe_NFE_AF","gnomADe_OTH_AF",
                   "gnomADe_SAS_AF","gnomADg_AF","gnomADg_AFR_AF","gnomADg_AMI_AF","gnomADg_AMR_AF",
                   "gnomADg_ASJ_AF","gnomADg_EAS_AF","gnomADg_FIN_AF","gnomADg_MID_AF","gnomADg_NFE_AF",
                   "gnomADg_OTH_AF","gnomADg_SAS_AF","MAX_AF","MAX_AF_POPS","CLIN_SIG","SOMATIC",
                   "PHENO","PUBMED","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE",
                   "TRANSCRIPTION_FACTORS","LoF","LoF_filter","LoF_flags","LoF_info"]

# NEW
def get_gnomAD_AF(csq, col_num):
    if type(csq)==float:
        return 0
    csqs = []
    for ind_csq in csq:
        af = ind_csq.split('|')[col_num]
        if af != '':
            csqs.append(af)
    csqs = list(set(csqs))
    if len(csqs)==0:
        return 0
    return csqs[0]

ultra_rare_vars_df['CSQ'] = ultra_rare_vars_df.CSQ.replace({'.':np.nan, None: np.nan})#.str.split(',')

try:
    n_csq_fields = len(ultra_rare_vars_df[~ultra_rare_vars_df.CSQ.isna()].CSQ.iloc[0][0].split('|'))

    if n_csq_fields==len(csq_columns_more):
        gnomad_af_str = 'gnomADe_AF'
        csq_columns = csq_columns_more
    elif n_csq_fields==len(csq_columns_less):
        gnomad_af_str = 'gnomAD_AF'
        csq_columns = csq_columns_less
    else:
        warnings.simplefilter("error")
        warnings.warn("CSQ fields are messed up!")

    ultra_rare_vars_df[gnomad_af_str] = ultra_rare_vars_df.CSQ.apply(get_gnomAD_AF, col_num=csq_columns.index(gnomad_af_str)).astype(float)
    ultra_rare_vars_df = ultra_rare_vars_df[ultra_rare_vars_df[gnomad_af_str]<=csq_af_threshold]
    cols_to_keep.append(gnomad_af_str)

except Exception as e:
    print(str(e))
    # pass

ultra_rare_vars_df['VarKey'] = ultra_rare_vars_df[['ID', 'REF', 'ALT']].astype(str).agg(':'.join, axis=1)
ultra_rare_vars_df[cols_to_keep].to_csv(f"{cohort_prefix}_ultra_rare_variants.tsv.gz", sep='\t', index=False)
