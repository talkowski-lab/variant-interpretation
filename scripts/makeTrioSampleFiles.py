import pandas as pd
import os
import sys
from google.cloud import storage

ped_uri = sys.argv[1]
cohort_prefix = sys.argv[2]
bucket_id = sys.argv[3]

storage_client = storage.Client(project='talkowski-sv-gnomad')
# bucket = storage_client.get_bucket(ped_uri.split('/')[2])
bucket = storage_client.get_bucket(bucket_id.replace('gs://', ''))

ped = pd.read_csv(ped_uri, sep='\t')

try:
    float(ped.columns[-1])  # no header
    ped = pd.read_csv(ped_uri, sep='\t', header=None)
    bucket.blob(f"resources/pedigrees/{cohort_prefix}_no_header.ped").upload_from_string(ped.to_csv(sep='\t', index=False, header=False), 'text/csv')
    ped.columns = ['FamID', 'IndividualID', 'FatherID', 'MotherID', 'Gender', 'Affected']
    bucket.blob(f"resources/pedigrees/{cohort_prefix}.ped").upload_from_string(ped.to_csv(sep='\t', index=False), 'text/csv')
except Exception as e:
    ped.columns = ['FamID', 'IndividualID', 'FatherID', 'MotherID', 'Gender', 'Affected']
    bucket.blob(f"resources/pedigrees/{cohort_prefix}.ped").upload_from_string(ped.to_csv(sep='\t', index=False), 'text/csv')
    # save ped no header
    bucket.blob(f"resources/pedigrees/{cohort_prefix}_no_header.ped").upload_from_string(ped.to_csv(sep='\t', index=False, header=False), 'text/csv')

trio = ped[(ped.FatherID != '0') & (ped.MotherID != '0')].iloc[:, :4]
trio['TrioID'] = trio['FamID'].astype(str) + '-' + trio['IndividualID'].astype(str)
trio.rename(columns={'IndividualID': 'SampleID'}, inplace=True)
trio = trio[["FamID","TrioID","SampleID","FatherID","MotherID"]]
bucket.blob(f"resources/metadata/{cohort_prefix}_trio_list.txt").upload_from_string(trio.to_csv(sep='\t', index=False), 'text/csv')
print('Trio csv saved')

trio.columns.name = 'Role'
trio.index = trio.FamID

sample_data = trio[['SampleID', 'FatherID', 'MotherID']].stack().reset_index()
sample_data.columns = ['FamID', 'Role', 'SampleID']
sample_data = sample_data.replace({'SampleID': 'child', 'MotherID': 'mother', 'FatherID': 'father'})
bucket.blob(f"resources/metadata/{cohort_prefix}_sample_list.txt").upload_from_string(sample_data.to_csv(sep='\t', index=False), 'text/csv')
print('Metadata csv saved')
