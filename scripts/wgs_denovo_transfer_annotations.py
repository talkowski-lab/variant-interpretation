import pandas as pd
import numpy as np
import sys
import ast
import os


trio_denovo_tsv = sys.argv[1]
merged_split_trio_tsv = sys.argv[2]

trio_denovo_df = pd.read_csv(trio_denovo_tsv, sep='\t')
trio_denovo_df.columns = trio_denovo_df.columns.str.replace('#', '')

merged_split_trio_df = pd.read_csv(merged_split_trio_tsv, sep='\t')

merged_split_trio_df['VarKey'] = merged_split_trio_df.locus.astype(str) + ':' + merged_split_trio_df.alleles.apply(ast.literal_eval).apply(pd.Series).apply(':'.join, axis=1) + ':' + merged_split_trio_df.SAMPLE.astype(str)
merged_split_trio_df.index = merged_split_trio_df.VarKey

trio_denovo_df['VarKey'] = trio_denovo_df[['CHROM', 'POS', 'REF', 'ALT']].astype(str).apply(':'.join, axis=1) + ':' + trio_denovo_df.SAMPLE.astype(str)
trio_denovo_df.index = trio_denovo_df.VarKey

merged_split_trio_df['in_step04'] = merged_split_trio_df.index.isin(trio_denovo_df.index)

merged_split_trio_df.to_csv(f"{os.path.basename(trio_denovo_tsv).split('.tsv')[0]}_annot.tsv", sep='\t', index=False)
