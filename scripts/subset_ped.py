import pandas as pd
import numpy as np
import sys

samples = sys.argv[1]
ped_uri = sys.argv[2]

ped = pd.read_csv(ped_uri, sep="\t")
samples = pd.read_csv(samples, header=None)[0]

ped.index = ped.sample_id
ped = ped.loc[samples]

for sample in ped.index:
    mother = ped.loc[sample, 'maternal_id']
    father = ped.loc[sample, 'paternal_id']
    if (mother == '0') & (father == '0'):
        continue
    try:
        ped.loc[mother]
        ped.loc[father]
    except Exception as e:
        # incomplete family
        ped = ped.drop(sample)

new_ped_filename = ped_uri.split('.ped')[0]+'_subset.ped'
ped.to_csv(new_ped_filename, sep="\t", index=False)
print(f"ped subset saved to {new_ped_filename}")
