import pandas as pd
import numpy as np
import sys
import os

samples = sys.argv[1]
ped_uri = sys.argv[2]

ped = pd.read_csv(ped_uri, sep="\t", dtype={i: str for i in range(4)})#.iloc[:,:6]
# ped.columns = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
samples = pd.read_csv(samples, header=None, dtype=str)[0]

ped.index = ped.sample_id
ped = ped.loc[np.intersect1d(ped.index, samples)]

fam_sizes = ped.family_id.value_counts().to_dict()
fathers = ped[ped.paternal_id!='0'].set_index('sample_id').paternal_id.to_dict()
mothers = ped[ped.maternal_id!='0'].set_index('sample_id').maternal_id.to_dict()

def get_sample_role(row):
    if fam_sizes[row.family_id]==1:
        role = 'Singleton'
    elif (row.maternal_id=='0') & (row.paternal_id=='0'):
        if (row.sample_id in fathers.values()):
            role = 'Father'
        elif (row.sample_id in mothers.values()):
            role = 'Mother'
        else:
            role = 'Unknown'
    elif (row.maternal_id=='-9') & (row.paternal_id=='-9'):
        role = 'Unknown'
    else:
        role = 'Proband'
    return role

ped['role'] = ped.apply(get_sample_role, axis=1)

# check 
good_fams = []
for family in ped.family_id.unique():
    fam_df = ped[ped.family_id==family].copy()
    probands = fam_df[fam_df.role=='Proband'].sample_id.tolist()
    mother = fam_df[fam_df.role=='Mother'].sample_id.tolist()
    father = fam_df[fam_df.role=='Father'].sample_id.tolist()
    if len(probands)==0:
        print(f"Family {family} has no probands!")      
        continue
    if len(mother)==0:
        print(f"Family {family} has no mother!")     
        continue
    if len(father)==0:    
        print(f"Family {family} has no father!")   
        continue
    mother, father = mother[0], father[0]
    
    if fam_df.loc[mother, 'sex']!=2:
        print(f"Mother {mother} in family {family} has wrong sex!")
        continue
    if fam_df.loc[father, 'sex']!=1:
        print(f"Father {father} in family {family} has wrong sex!")
        continue
    good_fams.append(family)

ped = ped[ped.family_id.isin(good_fams)].iloc[:, :6]

new_ped_filename = os.path.basename(ped_uri).split('.ped')[0]+'_subset.ped'
ped.to_csv(new_ped_filename, sep="\t", index=False)
print(f"ped subset saved to {new_ped_filename}")
