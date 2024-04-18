import hail as hl
import pandas as pd
import numpy as np
import sys
import ast
import warnings
import os

build = 'GRCh38'
vcf_file = sys.argv[1]
relatedness_qc = sys.argv[2]
ped_sex_qc = sys.argv[3]

hl.init()

mt = hl.import_vcf(vcf_file, reference_genome='GRCh38', force_bgz=True, call_fields=[], array_elements_required=False)
rel_df = pd.read_csv(relatedness_qc, sep='\t').set_index('sample_id')
ped_qc = pd.read_csv(ped_sex_qc, sep='\t').set_index('sample_id')

sample_qc = hl.sample_qc(mt).cols().flatten().to_pandas()

sample_qc.index = sample_qc.s
sample_qc.columns = sample_qc.columns.str.replace('sample_qc.','')

ped_qc['sex_error'] = ped_qc['sex']!=ped_qc['ped_sex']

rel_df['maternal_error'] = (rel_df.role=='Proband')&(~rel_df.mother_status.isin(['parent-child','ambiguous',np.nan]))
rel_df['paternal_error'] = (rel_df.role=='Proband')&(~rel_df.father_status.isin(['parent-child','ambiguous',np.nan]))

rel_df['missing_parent'] = (((rel_df[['maternal_id','paternal_id']]=='0')
            | (rel_df[['maternal_id','paternal_id']].isna())).sum(axis=1)==1)

sample_qc[['maternal_error','paternal_error','missing_parent']] = rel_df[['maternal_error','paternal_error','missing_parent']]
sample_qc['parental_error'] = sample_qc[['maternal_error','paternal_error']].any(axis=1)
sample_qc['both_parental_error'] = sample_qc[['maternal_error','paternal_error']].all(axis=1)
sample_qc['sex_error'] = ped_qc['sex_error']

sample_qc['any_error'] = sample_qc[['parental_error','sex_error']].any(axis=1)

n_het_thr = sample_qc.groupby('any_error').n_het.describe().loc[False, 'max']
samples_to_keep = sample_qc[(sample_qc.n_het<=n_het_thr)].index.tolist()

mt = mt.filter_cols(hl.array(samples_to_keep).contains(mt.s))

# clean-up: remove AC = 0 loci
mt = hl.variant_qc(mt)
mt = mt.filter_rows(mt.variant_qc.AC[1] > 0, keep = True)
mt = mt.drop('variant_qc')

sample_qc.to_csv(f"{os.path.basename(mt).split('.vcf')[0]}.sample_qc.txt", sep='\t', index=False)
hl.export_vcf(mt, os.path.basename(mt).split('.vcf')[0]+'_removed_outliers.vcf.bgz')