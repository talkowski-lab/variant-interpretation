import pandas as pd
import numpy as np
import hail as hl
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
import ast
import datetime

rel_ht_uri = sys.argv[1]
cohort_prefix = sys.argv[2]
ped_uri = sys.argv[3]
cores = sys.argv[4]  # string
mem = int(np.floor(float(sys.argv[5])))

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

rel = hl.read_table(rel_ht_uri)

ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
ped.columns = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']

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
    else:
        role = 'Proband'
    return role

ped['role'] = ped.apply(get_sample_role, axis=1)

ped_rels = {'i':[], 'j': [], 'ped_relationship': [], 'family_id': []}

for fam in ped.family_id.unique():
    fam_df = ped[ped.family_id==fam].reset_index()
    for i, row_i in fam_df.iterrows():
        for j, row_j in fam_df.iterrows():
            related = False
            if i==j:
                continue
                
            if (((row_i.paternal_id == row_j.sample_id)) | ((row_i.sample_id == row_j.paternal_id))) |\
            (((row_i.maternal_id == row_j.sample_id)) | ((row_i.sample_id == row_j.maternal_id))):                
                ped_rels['ped_relationship'].append('parent-child')
                related = True
            if (row_i.role not in ['Mother', 'Father']) & (row_j.role not in ['Mother', 'Father']):
                if (row_i.paternal_id==row_j.paternal_id) & (row_i.maternal_id==row_j.maternal_id):
                    ped_rels['ped_relationship'].append('siblings')
                else:
                    ped_rels['ped_relationship'].append('related_other')      
                related = True

            if related:
                ped_rels['i'].append(row_i.sample_id)
                ped_rels['j'].append(row_j.sample_id)
                ped_rels['family_id'].append(fam)


ped_rels_df = pd.DataFrame(ped_rels)

ped_rels_ht = hl.Table.from_pandas(ped_rels_df)

ped_rels_ht_merged = ped_rels_ht.annotate(i=ped_rels_ht.j, j=ped_rels_ht.i).key_by('i', 'j').union(ped_rels_ht.key_by('i','j'))

rel_merged = rel.key_by()
rel_merged = rel_merged.annotate(i=rel_merged.j, j=rel_merged.i).key_by('i', 'j').union(rel.key_by('i','j'))

related_in_ped = rel_merged.annotate(ped_relationship=ped_rels_ht_merged[rel_merged.key].ped_relationship)
related_in_ped = related_in_ped.filter(hl.is_missing(related_in_ped.ped_relationship), keep=False)

unrelated_in_ped = rel.anti_join(related_in_ped).annotate(ped_relationship='unrelated')

p = 0.05
only_related = unrelated_in_ped.filter(unrelated_in_ped.relationship!='unrelated')
downsampled_unrelated = unrelated_in_ped.filter(unrelated_in_ped.relationship=='unrelated').sample(p)

rel_total = related_in_ped.union(only_related).union(downsampled_unrelated)

rel_df = rel_total.to_pandas()

fig, ax = plt.subplots(1, 2, figsize=(12, 5));
fig.suptitle(cohort_prefix);

sns.scatterplot(rel_df, x='ibd0', y='kin', hue='relationship', s=16, ax=ax[0],
               hue_order=['parent-child', 'siblings', 'second degree relatives', 'duplicate/twins', 'ambiguous', 'unrelated'], 
               palette={'parent-child': 'mediumpurple', 'siblings': 'mediumseagreen', 'second degree relatives': 'skyblue', 
                        'duplicate/twins': 'indianred', 'ambiguous': 'sandybrown', 'unrelated': 'silver'});
ax[0].set_title("Inferred relationship from VCF");

sns.scatterplot(rel_df, x='ibd0', y='kin', hue='ped_relationship', s=16, ax=ax[1],
               hue_order=['parent-child', 'siblings', 'related_other', 'unrelated'],
               palette={'parent-child': 'mediumpurple', 'siblings': 'mediumseagreen', 'related_other': 'skyblue', 'unrelated': 'silver'});
ax[1].set_title("Relationship from pedigree");

plt.tight_layout();
plt.savefig(f"{cohort_prefix}_relatedness_ibd0_kinship.png");