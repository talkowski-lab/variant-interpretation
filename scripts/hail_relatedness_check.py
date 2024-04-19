import pandas as pd
import numpy as np
import hail as hl
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
import ast
import datetime
from gnomad.resources.grch38 import gnomad
from gnomad.sample_qc import relatedness

vcf_uri = sys.argv[1]
somalier_vcf = sys.argv[2]
cohort_prefix = sys.argv[3]
ped_uri = sys.argv[4]
cores = sys.argv[5]  # string
mem = int(np.floor(float(sys.argv[6])))
bucket_id = sys.argv[7]
score_table = sys.argv[8]

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

#split-multi
def split_multi_ssc(mt):
    mt = mt.annotate_rows(num_alleles = mt.alleles.size() ) # Add number of alleles at site before split
    # Now split
    sm = hl.split_multi(mt)
    pl = hl.or_missing(hl.is_defined(sm.PL),
                      (hl.range(0, 3).map(lambda i: hl.min(hl.range(0, hl.len(sm.PL))
       .filter(lambda j: hl.downcode(hl.unphased_diploid_gt_index_call(j), sm.a_index) == hl.unphased_diploid_gt_index_call(i))
       .map(lambda j: sm.PL[j])))))
    split_ds = sm.annotate_entries(GT = hl.downcode(sm.GT, sm.a_index),
                                   AD = hl.or_missing(hl.is_defined(sm.AD), [hl.sum(sm.AD) - sm.AD[sm.a_index], sm.AD[sm.a_index]]),
                                   PL = pl) 
        #GQ = hl.cond(hl.is_defined(pl[0]) & hl.is_defined(pl[1]) & hl.is_defined(pl[2]), hl.gq_from_pl(pl), sm.GQ) )
    mt = split_ds.drop('old_locus', 'old_alleles')
    return mt

mt = hl.import_vcf(vcf_uri, reference_genome='GRCh38', force_bgz=True, call_fields=[], array_elements_required=False)
mt = split_multi_ssc(mt)

# somalier sites
som_mt = hl.import_vcf(somalier_vcf, reference_genome='GRCh38', force_bgz=True, call_fields=[], array_elements_required=False)
mt = mt.semi_join_rows(som_mt.rows())

if score_table=='false':
    rel = hl.pc_relate(mt.GT, 0.01, k=10)
else:
    score_table_som = hl.read_table(score_table)
    mt = mt.annotate_cols(scores=score_table_som[mt.s].scores)
    rel = hl.pc_relate(mt.GT, 0.01, scores_expr=mt.scores)

rel = rel.annotate(relationship = relatedness.get_relationship_expr(rel.kin, rel.ibd0, rel.ibd1, rel.ibd2, 
                                                   first_degree_kin_thresholds=(0.19, 0.4), second_degree_min_kin=0.1, 
                                                   ibd0_0_max=0.2, ibd0_25_thresholds=(0.2, 0.425), ibd1_0_thresholds=(-0.15, 0.1), 
                                                   ibd1_50_thresholds=(0.275, 0.75), ibd1_100_min=0.75, ibd2_0_max=0.125, 
                                                   ibd2_25_thresholds=(0.1, 0.5), ibd2_100_thresholds=(0.75, 1.25))
)
rel = rel.key_by()
rel = rel.annotate(i=rel.i.s, j=rel.j.s).key_by('i','j')

# filename = f"{bucket_id}/hail/relatedness/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))}/{cohort_prefix}_pc_relate.ht"
# rel.write(filename)
# pd.Series([filename]).to_csv('mt_uri.txt', header=None, index=False)

ped = pd.read_csv(ped_uri, sep='\t').iloc[:, :6]
ped.columns = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
ped.index = ped.sample_id
try:
    ped['sex'] = ped.sex.astype(str).str.lower().replace({'male': 1, 'female': 2})
    ped.loc[~ped.sex.isin([1,2,'1','2']), 'sex'] = 0
    ped['sex'] = ped.sex.astype(int)
except Exception as e:
    print(str(e))
    pass
ped[['family_id', 'sample_id', 'paternal_id', 'maternal_id']] = ped[['family_id', 'sample_id', 'paternal_id', 'maternal_id']].astype(str)

# subset ped to samples in VCF
vcf_samps = mt.s.collect()
ped = ped[ped.sample_id.isin(vcf_samps)].copy()

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

ped_ht = hl.Table.from_pandas(ped)

dad_temp = ped_ht.key_by('sample_id','paternal_id')
dad_temp2 = ped_ht.key_by('paternal_id','sample_id')

all_dad_ht = dad_temp.join(rel.semi_join(dad_temp)).key_by().union(dad_temp2.join(rel.semi_join(dad_temp2)).key_by(), unify=True)
all_dad_ht = all_dad_ht.annotate(father_status=all_dad_ht.relationship).drop('relationship')

mom_temp = ped_ht.key_by('sample_id','maternal_id')
mom_temp2 = ped_ht.key_by('maternal_id','sample_id')

all_mom_ht = mom_temp.join(rel.semi_join(mom_temp)).key_by().union(mom_temp2.join(rel.semi_join(mom_temp2)).key_by(), unify=True)
all_mom_ht = all_mom_ht.annotate(mother_status=all_mom_ht.relationship).drop('relationship')

mom_df = all_mom_ht.to_pandas()
dad_df = all_dad_ht.to_pandas()

rename_cols = ['kin', 'ibd0', 'ibd1', 'ibd2']
mom_df = mom_df.rename({col: 'mother_'+col for col in rename_cols}, axis=1).copy()
dad_df = dad_df.rename({col: 'father_'+col for col in rename_cols}, axis=1).copy()

all_df = mom_df.merge(dad_df, how='outer')

# get duplicates
ped['sample_rank'] = 4 - (ped.sex.isin([1,2]).astype(int) + (ped.paternal_id!='0').astype(int) + (ped.maternal_id!='0').astype(int))

for s in np.setdiff1d(vcf_samps, ped.sample_id):
    ped.at[s, 'sample_id'] = s
    ped.loc[s, 'sample_rank'] = 4

# smaller rank is better
rank_ht = hl.Table.from_pandas(ped[['sample_id','sample_rank']]).key_by('sample_id')

merged_rel_df = pd.concat([ped, all_df.set_index('sample_id')], axis=1)
merged_rel_df = merged_rel_df.loc[:,~merged_rel_df.columns.duplicated()].copy()

dups = relatedness.get_duplicated_samples(rel, i_col='i', j_col='j', rel_col='relationship')
if len(dups)>0:
    dup_ht = relatedness.get_duplicated_samples_ht(dups, rank_ht,  rank_ann='sample_rank')
    dup_df = dup_ht.to_pandas()
    merged_rel_df['duplicate_samples'] = merged_rel_df.sample_id.map(dup_df.set_index('kept').filtered.to_dict())
else:
    merged_rel_df['duplicate_samples'] = np.nan

merged_rel_df.to_csv(f"{cohort_prefix}_relatedness_qc.ped", sep='\t', index=False)

# annotate ped
ped_rels = {'i':[], 'j': [], 'ped_relationship': [], 'family_id': []}

for fam in ped.family_id.unique():
    fam_df = ped[ped.family_id==fam].reset_index(drop=True)
    for i, row_i in fam_df.iterrows():
        for j, row_j in fam_df.iterrows():
            if (i==j) | ((row_i.role in ['Mother','Father']) & (row_j.role in ['Mother','Father'])):
                continue

            if (((row_i.paternal_id == row_j.sample_id)) | ((row_i.sample_id == row_j.paternal_id))) |\
            (((row_i.maternal_id == row_j.sample_id)) | ((row_i.sample_id == row_j.maternal_id))):                
                ped_rels['ped_relationship'].append('parent-child')
            elif (row_i.paternal_id==row_j.paternal_id) & (row_i.maternal_id==row_j.maternal_id):
                    ped_rels['ped_relationship'].append('siblings')
            else:
                ped_rels['ped_relationship'].append('related_other')      

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
rel_df.to_csv(f"{cohort_prefix}_kinship.tsv", sep='\t', index=False)