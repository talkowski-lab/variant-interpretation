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

final_output = final_output[(final_output.PL_sample!='.')&(~final_output.PL_sample.str.contains('\.'))]
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

gnomad_af_str = 'gnomADe_AF'
final_output = final_output[final_output[gnomad_af_str].replace({np.nan: 0})<=csq_af_threshold]

final_output.to_csv(os.path.basename(vcf_metrics_uri).split('.tsv')[0]+'_filtered.tsv.gz', sep='\t', index=False)