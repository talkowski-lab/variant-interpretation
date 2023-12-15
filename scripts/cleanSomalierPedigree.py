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

somalier = somalier.replace({'-9': '0'})
new_ped = ped.copy()

new_ped['sex_error'] = somalier.sex!=ped.sex
new_ped['family_error'] = somalier.family_id!=ped.family_id
new_ped['maternal_error'] = somalier.maternal_id!=ped.maternal_id
new_ped['paternal_error'] = somalier.paternal_id!=ped.paternal_id

ped[~(new_ped.sex_error | new_ped.family_error | new_ped.maternal_error | new_ped.paternal_error)].reset_index(drop=True).sort_values('sample_id').to_csv(f"{cohort_prefix}_ped_corrected.ped", sep='\t', index=False)

errors_df = pd.DataFrame({'sex_error':new_ped[new_ped.sex_error].sample_id,
                        'family_error': new_ped[new_ped.family_error].sample_id,
                        'maternal_error':new_ped[new_ped.maternal_error].sample_id,
                        'paternal_error':new_ped[new_ped.paternal_error].sample_id}).reset_index(drop=True)

errors_df.to_csv(f"{cohort_prefix}_somalier_errors.tsv", index=False, sep="\t", na_rep="")