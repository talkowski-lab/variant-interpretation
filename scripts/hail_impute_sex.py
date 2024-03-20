import pandas as pd
import numpy as np
import hail as hl
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
import ast

vcf_uri = sys.argv[1]
bed_file = sys.argv[2]
cohort_prefix = sys.argv[3]
ped_uri = sys.argv[4]
cores = sys.argv[5]  # string
mem = int(np.floor(float(sys.argv[6])))

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

mt = hl.import_vcf(vcf_uri, reference_genome='GRCh38', force_bgz=True, array_elements_required=False)
mt = split_multi_ssc(mt)

# somalier sites
intervals = hl.import_bed(bed_file, reference_genome='GRCh38')
mt = mt.filter_rows(hl.is_defined(intervals[mt.locus]))

imputed_sex = hl.impute_sex(mt.GT)
sex_df = imputed_sex.to_pandas()

sex_df['sex'] = sex_df.is_female.map({False: 1, True: 2})
sex_df.loc[sex_df.is_female.isna(), 'sex'] = 0

sex_df.index = sex_df.s.astype(str)
sex_df.index.name = 'sample_id'

all_qc = hl.sample_qc(mt)
all_qc_df = all_qc.cols().flatten().to_pandas()
all_qc_df = all_qc_df.rename({'s': 'sample_id'}, axis=1).set_index('sample_id')

gt_qc = hl.sample_qc(mt.filter_entries(hl.is_missing(mt.GT), keep=False))
gt_qc_df = gt_qc.cols().flatten().to_pandas()
gt_qc_df = gt_qc_df.rename({'s': 'sample_id'}, axis=1).set_index('sample_id')

sex_chrX_qc = hl.sample_qc(mt.filter_rows(mt.locus.in_x_nonpar()))
sex_chrX_qc_df = sex_chrX_qc.cols().flatten().to_pandas()
sex_chrX_qc_df = sex_chrX_qc_df.rename({'s': 'sample_id'}, axis=1).set_index('sample_id')

sex_chrY_qc = hl.sample_qc(mt.filter_rows(mt.locus.in_y_nonpar()))
sex_chrY_qc_df = sex_chrY_qc.cols().flatten().to_pandas()
sex_chrY_qc_df = sex_chrY_qc_df.rename({'s': 'sample_id'}, axis=1).set_index('sample_id')

sex_chrY_qc_df.columns = sex_chrY_qc_df.columns.str.replace('sample_qc', 'chrY')
sex_chrX_qc_df.columns = sex_chrX_qc_df.columns.str.replace('sample_qc', 'chrX')
all_qc_df.columns = all_qc_df.columns.str.replace('sample_qc', 'all')
gt_qc_df.columns = gt_qc_df.columns.str.replace('sample_qc', 'gt')

sample_qc_df = pd.concat([sex_chrX_qc_df, sex_chrY_qc_df, all_qc_df, gt_qc_df], axis=1)
sample_qc_df['hail_sex'] = sex_df.sex.astype('category')

sample_qc_df['Y_ploidy'] = 2 * sample_qc_df['chrY.dp_stats.mean'] / sample_qc_df['gt.dp_stats.mean']
sample_qc_df['X_ploidy'] = 2 * sample_qc_df['chrX.dp_stats.mean'] / sample_qc_df['gt.dp_stats.mean']

# sample_qc_df = sample_qc_df.loc[:,~sample_qc_df.columns.duplicated()].copy()

fig, ax = plt.subplots(2, 2, figsize=(12, 10));
sns.scatterplot(data=sample_qc_df, x='chrX.n_hom_var', y='chrX.n_het', hue='hail_sex', ax=ax[0][0]);
sns.scatterplot(data=sample_qc_df, x='chrX.n_hom_var', y='chrY.dp_stats.mean', hue='hail_sex', ax=ax[0][1]);
sns.scatterplot(data=sample_qc_df, x='chrX.dp_stats.mean', y='chrY.dp_stats.mean', hue='hail_sex', ax=ax[1][0]);
sns.scatterplot(data=sample_qc_df, x='X_ploidy', y='Y_ploidy', hue='hail_sex', ax=ax[1][1]);
plt.savefig(f"{cohort_prefix}_sex_qc_imputed.png");

ped = pd.read_csv(ped_uri, sep='\t').iloc[:, :6]
base_cols = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
ped.columns = base_cols
ped = ped.set_index('sample_id')

ped_qc = pd.concat([ped, sample_qc_df], axis=1).reset_index()[base_cols + np.setdiff1d(sample_qc_df.columns, base_cols).tolist()]
ped_qc.to_csv(f"{cohort_prefix}_sex_qc.ped", sep='\t', index=False)

fig, ax = plt.subplots(2, 2, figsize=(12, 10));
sns.scatterplot(data=ped_qc, x='chrX.n_hom_var', y='chrX.n_het', hue='sex', ax=ax[0][0]);
sns.scatterplot(data=ped_qc, x='chrX.n_hom_var', y='chrY.dp_stats.mean', hue='sex', ax=ax[0][1]);
sns.scatterplot(data=ped_qc, x='chrX.dp_stats.mean', y='chrY.dp_stats.mean', hue='sex', ax=ax[1][0]);
sns.scatterplot(data=ped_qc, x='X_ploidy', y='Y_ploidy', hue='sex', ax=ax[1][1]);
plt.savefig(f"{cohort_prefix}_sex_qc_ped.png");
