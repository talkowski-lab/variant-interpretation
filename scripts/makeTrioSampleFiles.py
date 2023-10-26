import pandas as pd
import os
import sys

ped_uri = sys.argv[1].split('/')[-1]
cohort_prefix = sys.argv[2]

ped = pd.read_csv(ped_uri, sep='\t')
trio = ped[(ped.FatherID != '0') & (ped.MotherID != '0')].iloc[:, :4]
trio['TrioID'] = trio['FamID'] + '-' + trio['IndividualID']
trio.rename(columns={'IndividualID': 'SampleID'}, inplace=True)
trio.to_csv(f"{cohort_prefix}_trio_list.txt", sep='\t', index=False)
print('Trio csv saved')

trio.columns.name = 'Role'
trio.index = trio.FamID

sample_data = trio.iloc[:, 1:-1].stack().reset_index()
sample_data.columns = ['FamID', 'Role', 'SampleID']
sample_data = sample_data.replace({'SampleID': 'child', 'MotherID': 'mother', 'FatherID': 'father'})
sample_data.to_csv(f"{cohort_prefix}_sample_list.txt", sep='\t', index=False)
print('Metadata csv saved')