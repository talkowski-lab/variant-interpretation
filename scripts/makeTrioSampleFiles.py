import pandas as pd
import os
import sys
from google.cloud import storage

ped_uri = sys.argv[1]
cohort_prefix = sys.argv[2]
print(ped_uri)

storage_client = storage.Client(project='talkowski-sv-gnomad')
bucket = storage_client.get_bucket(ped_uri.split('gs://')[1].split('/')[0])

ped = pd.read_csv(ped_uri, sep='\t')
trio = ped[(ped.FatherID != '0') & (ped.MotherID != '0')].iloc[:, :4]
trio['TrioID'] = trio['FamID'] + '-' + trio['IndividualID']
trio.rename(columns={'IndividualID': 'SampleID'}, inplace=True)
# trio.to_csv(f"{cohort_prefix}_trio_list.txt", sep='\t', index=False)
bucket.blob(f"{cohort_prefix}_trio_list.txt").upload_from_string(trio.to_csv(sep='\t', index=False), 'text/csv')
print('Trio csv saved')

trio.columns.name = 'Role'
trio.index = trio.FamID

sample_data = trio.iloc[:, 1:-1].stack().reset_index()
sample_data.columns = ['FamID', 'Role', 'SampleID']
sample_data = sample_data.replace({'SampleID': 'child', 'MotherID': 'mother', 'FatherID': 'father'})
# sample_data.to_csv(f"{cohort_prefix}_sample_list.txt", sep='\t', index=False)
bucket.blob(f"{cohort_prefix}_sample_list.txt").upload_from_string(sample_data.to_csv(sep='\t', index=False), 'text/csv')
print('Metadata csv saved')

trio_uri = f"{ped_uri.split('/')[0]}/{cohort_prefix}_trio_list.txt"
meta_uri = f"{ped_uri.split('/')[0]}/{cohort_prefix}_sample_list.txt"
return (trio_uri, meta_uri)
