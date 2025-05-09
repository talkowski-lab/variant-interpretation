import pandas as pd
import numpy as np
import hail as hl
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
import ast

vcf_uri = sys.argv[1]
somalier_vcf = sys.argv[2]
cohort_prefix = sys.argv[3]
ped_uri = sys.argv[4]
cores = sys.argv[5]  # string
mem = int(np.floor(float(sys.argv[6])))
genome_build = sys.argv[7]
split_multi = ast.literal_eval(sys.argv[8].capitalize())

hl.init(min_block_size=128, 
        local=f"local[*]", 
        spark_conf={
                    "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.speculation": 'true'
                    }, 
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

#split-multi
def split_multi_ssc(mt):
    mt = mt.annotate_rows(num_alleles = mt.alleles.size() ) # Add number of alleles at site before split
    # only split variants that aren't already split
    bi = mt.filter_rows(hl.len(mt.alleles) == 2)
    bi = bi.annotate_rows(a_index=1, was_split=False, old_locus=bi.locus, old_alleles=bi.alleles)
    multi = mt.filter_rows(hl.len(mt.alleles) > 2)
    # Now split
    split = hl.split_multi(multi, permit_shuffle=True)
    sm = split.union_rows(bi)
    # sm = hl.split_multi(mt, permit_shuffle=True)
    if 'PL' in list(mt.entry.keys()):
        pl = hl.or_missing(hl.is_defined(sm.PL),
                        (hl.range(0, 3).map(lambda i: hl.min(hl.range(0, hl.len(sm.PL))
        .filter(lambda j: hl.downcode(hl.unphased_diploid_gt_index_call(j), sm.a_index) == hl.unphased_diploid_gt_index_call(i))
        .map(lambda j: sm.PL[j])))))
        sm = sm.annotate_entries(PL = pl)
    split_ds = sm.annotate_entries(GT = hl.downcode(sm.GT, sm.a_index),
                                   AD = hl.or_missing(hl.is_defined(sm.AD), [hl.sum(sm.AD) - sm.AD[sm.a_index], sm.AD[sm.a_index]])
                                   ) 
        #GQ = hl.cond(hl.is_defined(pl[0]) & hl.is_defined(pl[1]) & hl.is_defined(pl[2]), hl.gq_from_pl(pl), sm.GQ) )
    mt = split_ds.drop('old_locus', 'old_alleles')
    return mt

mt = hl.import_vcf(vcf_uri, reference_genome=genome_build, force_bgz=True, call_fields=[], array_elements_required=False)
if split_multi:
    mt = split_multi_ssc(mt)

# somalier sites
som_mt = hl.import_vcf(somalier_vcf, reference_genome=genome_build, force_bgz=True, call_fields=[], array_elements_required=False)
mt = mt.semi_join_rows(som_mt.rows())

mt = mt.annotate_entries(AB=mt.AD[1]/hl.sum(mt.AD))
mt = mt.annotate_entries(GT=hl.or_missing(mt.DP!=0, mt.GT))
# only consider hets with 0.3 <= AB <= 0.7
mt = mt.annotate_entries(GT=hl.if_else(mt.GT.is_het(),
                         hl.or_missing((mt.AB>=0.3) & (mt.AB<=0.7), mt.GT), mt.GT))

all_qc = hl.sample_qc(mt)
all_qc_df = all_qc.cols().flatten().to_pandas()
all_qc_df = all_qc_df.rename({'s': 'sample_id'}, axis=1).set_index('sample_id')

gt_mt = mt.filter_entries(hl.is_missing(mt.GT), keep=False)
gt_mt = gt_mt.filter_entries(gt_mt.DP>=7)
gt_qc = hl.sample_qc(gt_mt)
gt_qc_df = gt_qc.cols().flatten().to_pandas()
gt_qc_df = gt_qc_df.rename({'s': 'sample_id'}, axis=1).set_index('sample_id')

x_mt = mt.filter_rows(mt.locus.in_x_nonpar())
x_mt = x_mt.filter_entries(x_mt.DP>=7)
sex_chrX_qc = hl.sample_qc(x_mt)
sex_chrX_qc_df = sex_chrX_qc.cols().flatten().to_pandas()
sex_chrX_qc_df = sex_chrX_qc_df.rename({'s': 'sample_id'}, axis=1).set_index('sample_id')

y_mt = mt.filter_rows(mt.locus.in_y_nonpar())
sex_chrY_qc = hl.sample_qc(y_mt)
sex_chrY_qc_df = sex_chrY_qc.cols().flatten().to_pandas()
sex_chrY_qc_df = sex_chrY_qc_df.rename({'s': 'sample_id'}, axis=1).set_index('sample_id')

sex_chrY_qc_df.columns = sex_chrY_qc_df.columns.str.replace('sample_qc', 'chrY')
sex_chrX_qc_df.columns = sex_chrX_qc_df.columns.str.replace('sample_qc', 'chrX')
all_qc_df.columns = all_qc_df.columns.str.replace('sample_qc', 'all')
gt_qc_df.columns = gt_qc_df.columns.str.replace('sample_qc', 'gt')

sample_qc_df = pd.concat([sex_chrX_qc_df, sex_chrY_qc_df, all_qc_df, gt_qc_df], axis=1)

sample_qc_df['Y_ploidy'] = 2 * sample_qc_df['chrY.dp_stats.mean'] / sample_qc_df['gt.dp_stats.mean']
sample_qc_df['X_ploidy'] = 2 * sample_qc_df['chrX.dp_stats.mean'] / sample_qc_df['gt.dp_stats.mean']

# Somalier's sex prediction criteria
def predict_sex(row):
    sex = -9
    if row['chrX.n_hom_var'] > 0:
        if (row['chrX.n_het'] / row['chrX.n_hom_var'] < 0.05) & (row['chrX.n_called'] > 10):
            sex = 1
        elif (row['chrX.n_het'] / row['chrX.n_hom_var'] > 0.4) & (row['chrX.n_called'] > 10):
            sex = 2
    if (row['chrY.n_called']>0) & (row.sex==1) & (row['Y_ploidy'] < 0.4):
        sex = -1  # apparent loss of Y with low het-ratio on X chromosome
    if (row['chrY.n_called']>0) & (row.sex==2) & (row['Y_ploidy'] > 0.4):
        sex = -2  # apparent Y with high het-ratio on X chromosome
    return sex

ped = pd.read_csv(ped_uri, sep='\t').iloc[:, :6]
base_cols = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
ped.columns = base_cols
ped = ped.set_index('sample_id')

# subset ped to samples in VCF
vcf_samps = mt.s.collect()
ped = ped[ped.index.isin(vcf_samps)].copy()

ped_qc = pd.concat([ped, sample_qc_df], axis=1).reset_index()
ped_qc = ped_qc.fillna(np.nan)
ped_qc['ped_sex'] = ped_qc.sex
ped_qc['sex'] = ped_qc.apply(predict_sex, axis=1).astype('category')
ped_qc = ped_qc[base_cols + ['ped_sex'] + np.setdiff1d(sample_qc_df.columns, base_cols).tolist()]

for col in ped_qc.columns.tolist()[:6] + ['ped_sex']:
    ped_qc[col] = ped_qc[col].replace({np.nan: -9})
ped_qc.to_csv(f"{cohort_prefix}_sex_qc.ped", sep='\t', index=False)

sex_str = []
for sex_type in ped_qc.ped_sex.unique():
    try:
        float(sex_type)
    except:
        sex_str.append(sex_type)

ped_qc['ped_sex'] = ped_qc.ped_sex.replace({np.nan: -9} | {sex_type: -9 for sex_type in sex_str}).astype(int).astype('category')
ped_qc['sex'] = ped_qc.sex.astype('category')

fig = plt.figure(layout="constrained", figsize=(22, 10))
subfigs = fig.subfigures(1, 2, wspace=0.07, width_ratios=[1, 1])
ax = subfigs[0].subplots(2, 2);
subfigs[0].suptitle(f"{cohort_prefix}, Inferred Sex");
sns.scatterplot(data=ped_qc, x='chrX.n_hom_var', y='chrX.n_het', hue='sex',
               hue_order=np.intersect1d(ped_qc.sex.cat.categories, [2, 1, -1, -2, -9]), 
               palette={2: 'mediumpurple', 1: 'mediumseagreen', -1: 'skyblue', 
                        -2: 'palevioletred', -9: 'slategray'}, ax=ax[0][0]);
sns.scatterplot(data=ped_qc, x='chrX.n_hom_var', y='chrY.dp_stats.mean', hue='sex', 
               hue_order=np.intersect1d(ped_qc.sex.cat.categories, [2, 1, -1, -2, -9]), 
               palette={2: 'mediumpurple', 1: 'mediumseagreen', -1: 'skyblue', 
                        -2: 'palevioletred', -9: 'slategray'}, ax=ax[0][1]);
sns.scatterplot(data=ped_qc, x='chrX.dp_stats.mean', y='chrY.dp_stats.mean', hue='sex', 
               hue_order=np.intersect1d(ped_qc.sex.cat.categories, [2, 1, -1, -2, -9]), 
               palette={2: 'mediumpurple', 1: 'mediumseagreen', -1: 'skyblue', 
                        -2: 'palevioletred', -9: 'slategray'}, ax=ax[1][0]);
sns.scatterplot(data=ped_qc, x='X_ploidy', y='Y_ploidy', hue='sex', 
               hue_order=np.intersect1d(ped_qc.sex.cat.categories, [2, 1, -1, -2, -9]), 
               palette={2: 'mediumpurple', 1: 'mediumseagreen', -1: 'skyblue', 
                        -2: 'palevioletred', -9: 'slategray'}, ax=ax[1][1]);

ax = subfigs[1].subplots(2, 2);
subfigs[1].suptitle(f"{cohort_prefix}, Pedigree Sex");
sns.scatterplot(data=ped_qc, x='chrX.n_hom_var', y='chrX.n_het', hue='ped_sex', 
               hue_order=np.intersect1d(ped_qc.ped_sex.cat.categories, [2, 1, -1, -2, -9]), 
               palette={2: 'mediumpurple', 1: 'mediumseagreen', -1: 'skyblue', 
                        -2: 'palevioletred', -9: 'slategray'}, ax=ax[0][0]);
sns.scatterplot(data=ped_qc, x='chrX.n_hom_var', y='chrY.dp_stats.mean', hue='ped_sex', 
               hue_order=np.intersect1d(ped_qc.ped_sex.cat.categories, [2, 1, -1, -2, -9]), 
               palette={2: 'mediumpurple', 1: 'mediumseagreen', -1: 'skyblue', 
                        -2: 'palevioletred', -9: 'slategray'}, ax=ax[0][1]);
sns.scatterplot(data=ped_qc, x='chrX.dp_stats.mean', y='chrY.dp_stats.mean', hue='ped_sex', 
               hue_order=np.intersect1d(ped_qc.ped_sex.cat.categories, [2, 1, -1, -2, -9]), 
               palette={2: 'mediumpurple', 1: 'mediumseagreen', -1: 'skyblue', 
                        -2: 'palevioletred', -9: 'slategray'}, ax=ax[1][0]);
sns.scatterplot(data=ped_qc, x='X_ploidy', y='Y_ploidy', hue='ped_sex', 
               hue_order=np.intersect1d(ped_qc.ped_sex.cat.categories, [2, 1, -1, -2, -9]), 
               palette={2: 'mediumpurple', 1: 'mediumseagreen', -1: 'skyblue', 
                        -2: 'palevioletred', -9: 'slategray'}, ax=ax[1][1]);
plt.savefig(f"{cohort_prefix}_sex_qc.png");