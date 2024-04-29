import datetime
import pandas as pd
import hail as hl
import numpy as np
import pandas as pd
import ast
import sys
import os

de_novo_merged = sys.argv[1]
cohort_prefix = sys.argv[2]
ped_uri = sys.argv[3]
loeuf_file = sys.argv[4]
cores = sys.argv[5]
mem = int(np.floor(float(sys.argv[6])))
bucket_id = sys.argv[7]

prefix = os.path.basename(filtered_mt).split('_wes_denovo_basic_filtering.mt')[0]

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

df = pd.read_csv(de_novo_merged, sep='\t')
df['VarKey'] = df[['ID', 'proband.s']].astype(str).agg(':'.join, axis=1)
df['Consequence'] = df.Consequence.replace({np.nan: '[]'}).apply(ast.literal_eval).agg(','.join)

# Filter to only hets
df = df[~df['proband_entry.GT'].isin(['1/1','1|1'])]

# Filter out LOW confidence
df = df[df.confidence!='LOW']

# Filter SNPs on VQSLOD
df = df[(df.isIndel)|(df.VQSLOD >= -20)]   # TODO: input

# Filter indels on VQSLOD
df = df[(df.isSNV)|(df.VQSLOD >= -2)]   # TODO: input

# Set frequency threshold
MAF_thresh = 0.001  # TODO: input

# Filter on dataset AF
df = df[df.cohort_AF<= MAF_thresh]

# Filter on gnomAD
df = df[(df.gnomad_non_neuro_AF.isna())|(df.gnomad_non_neuro_AF<=MAF_thresh)]

# Impose allele balance filter
df = df[df['proband_entry.AB']>=0.3]  # TODO
df = df[(df['mother_entry.AB']<=0.05)&(df['father_entry.AB']<=0.05)] 

# Pick one variant per gene per sample
df['SAMPLE_GENE'] = df[['proband.s', 'SYMBOL']].astype(str).agg('-'.join, axis=1)
df = df.sort_values(['SAMPLE_GENE','csq_score'], ascending=False).groupby('SAMPLE_GENE').head(1).reset_index(drop=True)