import pandas as pd
import os
import sys
from google.cloud import storage

ped_uri = sys.argv[1]
cohort_prefix = sys.argv[2]

storage_client = storage.Client(project='talkowski-sv-gnomad')
bucket = storage_client.get_bucket(ped_uri.split('/')[2])

ped = pd.read_csv(ped_uri, sep='\t')
try:
    ped['IndividualID'] = ped['IndividualID'].str.split('__').str[1].str.upper()
    ped['FatherID'] = ped['FatherID'].str.split('__').str[1].str.upper().fillna(0)
    ped['MotherID'] = ped['MotherID'].str.split('__').str[1].str.upper().fillna(0)
    bucket.blob(f"resources/pedigrees/{ped_uri.split('/')[-1]}").upload_from_string(ped.to_csv(sep='\t', index=False), 'text/csv')
    print('Updated pedigree csv saved')
except Exception as e:
    pass

trio = ped[(ped.FatherID != '0') & (ped.MotherID != '0')].iloc[:, :4]
trio['TrioID'] = trio['FamID'] + '-' + trio['IndividualID']
trio.rename(columns={'IndividualID': 'SampleID'}, inplace=True)
bucket.blob(f"resources/metadata/{cohort_prefix}_trio_list.txt").upload_from_string(trio.to_csv(sep='\t', index=False), 'text/csv')
print('Trio csv saved')

trio.columns.name = 'Role'
trio.index = trio.FamID

sample_data = trio.iloc[:, 1:-1].stack().reset_index()
sample_data.columns = ['FamID', 'Role', 'SampleID']
sample_data = sample_data.replace({'SampleID': 'child', 'MotherID': 'mother', 'FatherID': 'father'})
bucket.blob(f"resources/metadata/{cohort_prefix}_sample_list.txt").upload_from_string(sample_data.to_csv(sep='\t', index=False), 'text/csv')
print('Metadata csv saved')
