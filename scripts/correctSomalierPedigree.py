import pandas as pd
import numpy as np
import os
import sys

samples_uri = sys.argv[1]
ped_uri = sys.argv[2]
cohort_prefix = sys.argv[3]
subset_ped = bool(sys.argv[4].upper())

ped = pd.read_csv(ped_uri, sep='\t')
try:
    float(ped.columns[-1])  # no header
    ped = pd.read_csv(ped_uri, sep='\t', header=None)
except Exception as e:
    pass

ped.columns = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
ped.index = ped.sample_id

somalier = pd.read_csv(samples_uri, sep='\t')
somalier.index = somalier.sample_id

if subset_ped:
    ped = ped.loc[np.intersect1d(ped.index, somalier.index)]

somalier = somalier.loc[ped.index]
somalier.columns = somalier.columns.str.replace("#", "")

discrepant_samples = somalier.loc[(ped.iloc[:,:6] != somalier.iloc[:,:6]).any(axis=1)].sample_id.tolist()

parent_ids = ['paternal_id', 'maternal_id']
for samp in discrepant_samples:
    for i, parent_id in enumerate(parent_ids):
        other_parent = parent_ids[i-1]
        if somalier.loc[samp, parent_id] == '-9':
            print(f"{samp} {parent_id} changed")
            somalier.loc[samp, parent_id] = 0
            if somalier.loc[samp, other_parent] in somalier.sample_id.tolist(): 
                if samp in somalier.loc[somalier.loc[samp, other_parent], parent_ids].tolist():
                    print(f"{samp} {other_parent} changed")
                    somalier.loc[samp, other_parent] = 0
                else:
                    # not a complete trio/something strange happening
                    print(f"{samp} is strange, dropping")
                    somalier = somalier.drop(samp)
                    break

somalier.reset_index(drop=True).sort_values('sample_id').iloc[:,:6].to_csv(f"{cohort_prefix}_ped_corrected.ped", sep='\t', index=False)
