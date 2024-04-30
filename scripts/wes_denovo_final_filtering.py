import datetime
import pandas as pd
import numpy as np
import pandas as pd
import ast
import sys
import os

de_novo_merged = sys.argv[1]
cohort_prefix = sys.argv[2]
vqslod_cutoff_snv = int(sys.argv[3])
vqslod_cutoff_indel = int(sys.argv[4])
MAF_thresh = float(sys.argv[5])
cores = sys.argv[6]
mem = int(np.floor(float(sys.argv[7])))

df = pd.read_csv(de_novo_merged, sep='\t')
df['VarKey'] = df[['ID', 'proband.s']].astype(str).agg(':'.join, axis=1)
df['Consequence'] = df.Consequence.replace({np.nan: '[]'}).apply(ast.literal_eval).agg(','.join)

# Filter to only hets
df = df[~df['proband_entry.GT'].isin(['1/1','1|1'])]

# Filter out LOW confidence
df = df[df.confidence!='LOW']

# Filter SNPs on VQSLOD
df = df[(df.isIndel)|(df.VQSLOD >= vqslod_cutoff_snv)] 

# Filter indels on VQSLOD
df = df[(df.isSNV)|(df.VQSLOD >= vqslod_cutoff_indel)] 

# Set frequency threshold

# Filter on dataset AF
df = df[df.cohort_AF<= MAF_thresh]

# Filter on gnomAD
df = df[(df.gnomad_non_neuro_AF.isna())|(df.gnomad_non_neuro_AF<=MAF_thresh)]

# Pick one variant per gene per sample
df['SAMPLE_GENE'] = df[['proband.s', 'SYMBOL']].astype(str).agg('-'.join, axis=1)
df = df.sort_values(['SAMPLE_GENE','csq_score'], ascending=False).groupby('SAMPLE_GENE').head(1).reset_index(drop=True)

df.to_csv(cohort_prefix+'_de_novo_filtered_final.tsv', sep='\t', index=False)