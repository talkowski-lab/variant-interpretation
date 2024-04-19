import os
import pandas as pd
import sys

ped_uri = sys.argv[1]
ped = pd.read_csv(ped_uri, sep='\t', dtype={i: str for i in range(4)}).iloc[:,:6]
ped.columns = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
ped.index = ped.sample_id

sample = sys.argv[2]

parents = ped.loc[sample].iloc[2:4].tolist()

ped.loc[parents+[sample]].to_csv(f"{sample}.ped", sep='\t', index=False, header=None)