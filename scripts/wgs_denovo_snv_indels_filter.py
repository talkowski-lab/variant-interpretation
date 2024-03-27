import os
import sys
import pandas as pd
import numpy as np
import warnings
import ast

vcf_metrics_uri = sys.argv[1]
AC_threshold = int(sys.argv[2])
AF_threshold = float(sys.argv[3])
csq_af_threshold = float(sys.argv[4])

final_output = pd.read_csv(vcf_metrics_uri, sep='\t')

final_output = final_output[final_output.PL_sample!='.']
final_output[['PL_sample_0.0', 'PL_sample_0.1', 'PL_sample_1.1']] = final_output.PL_sample.str.split(",", expand=True).astype(int)

for samp in ['sample', 'mother', 'father']:
    final_output[f"DPC_{samp}"] = final_output[f"AD_{samp}"].apply(ast.literal_eval).apply(sum)
    final_output[f"AB_{samp}"] = final_output[f"AD_{samp}"].apply(ast.literal_eval).str[1] / final_output[f"DPC_{samp}"]

final_output['GQ_parent'] = final_output[['GQ_mother', 'GQ_father']].min(axis=1)
final_output['AB_parent'] = final_output[['AB_mother', 'AB_father']].min(axis=1)
final_output['DPC_parent'] = final_output[['DPC_mother', 'DPC_father']].min(axis=1)
final_output['VAF_parent'] = final_output[['VAF_mother', 'VAF_father']].min(axis=1)

final_output.index = final_output.VarKey
final_output = final_output[final_output.POLYX <= 10]
# final_output = final_output[final_output['VQSLOD']!='.']
if 'cohort_AF' not in final_output.columns:
    final_output['cohort_AF'] = final_output.cohort_AC / final_output.AN
final_output = final_output[(final_output.cohort_AC<=AC_threshold) | (final_output.cohort_AF<=AF_threshold)]
final_output = final_output[final_output['PL_sample_0.1']==0]

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

# OLD
def get_csq_max_af(csq):
    if csq=='.':
        return ''
    csq_df = pd.DataFrame(csq.split(','))[0].str.split('|', expand=True)      
    try:
        csq_df.columns = csq_columns_more
    except:
        csq_df.columns = csq_columns_less            
    return csq_df.MAX_AF.max()

# final_output.loc[:,'MAX_AF'] = final_output.CSQ.apply(lambda csq: get_csq_max_af(csq)).replace({'': 0}).astype(float)
# final_output['MAX_AF'] = final_output.MAX_AF.replace({np.nan: 0})
# final_output = final_output[final_output.MAX_AF<=csq_af_threshold]

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

final_output['CSQ'] = final_output.CSQ.replace({'.':np.nan}).str.split(',')
n_csq_fields = len(final_output[~final_output.CSQ.isna()].CSQ.iloc[0][0].split('|'))

if n_csq_fields==len(csq_columns_more):
    gnomad_af_str = 'gnomADe_AF'
    csq_columns = csq_columns_more
elif n_csq_fields==len(csq_columns_less):
    gnomad_af_str = 'gnomAD_AF'
    csq_columns = csq_columns_less
else:
    warnings.simplefilter("error")
    warnings.warn("CSQ fields are messed up!")

final_output[gnomad_af_str] = final_output.CSQ.apply(get_gnomAD_AF, col_num=csq_columns.index(gnomad_af_str)).astype(float)
final_output = final_output[final_output[gnomad_af_str]<=csq_af_threshold]

final_output.to_csv(os.path.basename(vcf_metrics_uri).split('.tsv')[0]+'_filtered.tsv', sep='\t', index=False)