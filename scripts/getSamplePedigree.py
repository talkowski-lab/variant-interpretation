import os
import pandas as pd
import sys

ped_uri = sys.argv[1]
ped = pd.read_csv(ped_uri, sep='\t', header=None)
ped.columns = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
ped.index = ped.sample_id

vcf_file = sys.argv[2]

sample = os.path.basename(vcf_file).split('_trio_')[1].replace('.vcf', '')
parents = ped.loc[sample].iloc[2:4].tolist()

ped.loc[parents+[sample]].to_csv(f"{sample}.ped", sep='\t', index=False, header=None)