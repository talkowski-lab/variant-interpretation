# from wgs_denovo_snv_indels_qc.ipynb

import pandas as pd
import hail as hl
import os
import sys

cohort_tsv = sys.argv[1]
num_vars_step03_dir = sys.argv[2]
num_vars_step04_dir = sys.argv[3]

cohorts = pd.read_csv(cohort_tsv, sep="\t")
cohorts.index = cohorts['entity:cohort_id']

num_variants_per_step = pd.DataFrame()

# step01 output
for cohort, vcf_file in cohorts.merged_preprocessed_vcf_file.dropna().to_dict().items():
    mt = hl.import_vcf(vcf_file, force_bgz=True, call_fields=[], reference_genome='GRCh38')
    num_variants_per_step.loc[cohort, 'tot_step01'] = mt.count_rows()

# step03 output
for cohort in num_variants_per_step.index:
    num_vars_trio = pd.read_csv(f"{num_vars_step03_dir}/{cohort}_split_trio_num_vars.tsv", 
                                sep='\t', header=None)
    num_vars_trio.columns = ['uri', 'num_vars']
    num_vars_trio['trio'] = num_vars_trio.uri.apply(os.path.basename).str.split('.').str[0]
    num_vars_trio.index = num_vars_trio.trio
    trio_stats = num_vars_trio.num_vars.describe()
    if trio_stats['min']<0:
        print(f"{cohort} has {(num_vars_trio.num_vars<0).sum()} trios with a negative number of variants after step03...")
    num_variants_per_step.loc[cohort, 'num_trios_step03'] = num_vars_trio.shape[0]
    num_variants_per_step.loc[cohort, 'tot_step03'] = num_vars_trio.num_vars.sum()
    num_variants_per_step.loc[cohort, 'avg_step03'] = trio_stats['mean'].round(2).astype(str) + ' (+-' + trio_stats['std'].round(2).astype(str) + ')'

# step04 output
for cohort in num_variants_per_step.index:
    num_vars_trio = pd.read_csv(f"{num_vars_step04_dir}/{cohort}_trio_denovo_num_vars.tsv", 
                                sep='\t', header=None)
    num_vars_trio.columns = ['uri', 'num_vars']
    num_vars_trio['trio'] = num_vars_trio.uri.apply(os.path.basename).str.split('.').str[0]
    num_vars_trio.index = num_vars_trio.trio
    trio_stats = num_vars_trio.num_vars.describe()
    if trio_stats['min']<0:
        print(f"{cohort} has {(num_vars_trio.num_vars<0).sum()} trios with a negative number of variants after step04...")
    num_variants_per_step.loc[cohort, 'num_trios_step04'] = num_vars_trio.shape[0]
    num_variants_per_step.loc[cohort, 'tot_step04'] = num_vars_trio.num_vars.sum()
    num_variants_per_step.loc[cohort, 'avg_step04'] = trio_stats['mean'].round(2).astype(str) + ' (+-' + trio_stats['std'].round(2).astype(str) + ')'

# step05 output
for cohort, tsv_uri in cohorts.vcf_metrics_tsv.dropna().to_dict().items():
    final_output = pd.read_csv(tsv_uri, sep='\t')
    final_output.index = final_output.ID
    num_variants_per_step.loc[cohort, 'num_trios_step05'] = final_output.SAMPLE.unique().size
    num_variants_per_step.loc[cohort, 'tot_step05'] = final_output.ID.unique().size 
    samp_value_counts = final_output.SAMPLE.value_counts()
    num_variants_per_step.loc[cohort, 'avg_step05'] = samp_value_counts.mean().round(2).astype(str) + ' (+-' + samp_value_counts.std().round(2).astype(str) + ')'

num_variants_per_step = num_variants_per_step[['num_trios_step03', 'num_trios_step04', 'num_trios_step05',
                                               'tot_step01', 'tot_step03', 'tot_step04', 'tot_step05',
                                               'avg_step03', 'avg_step04', 'avg_step05']]

num_variants_per_step.index.name = 'cohort'
num_variants_per_step.to_csv(f"{os.path.basename(cohort_tsv).split('.tsv')[0]}_total_avg_num_vars_per_trio.tsv", sep='\t')