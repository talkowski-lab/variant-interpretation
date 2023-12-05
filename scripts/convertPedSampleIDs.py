import os
import sys
import pandas as pd
from google.cloud import storage

sample_tsv_uri = sys.argv[1]  # "gs://fc-71e715ea-2fb8-4a20-8560-15ed867dcc7d/resources/metadata/GMKF_collaborator_id_sample_id.tsv"
ped_uri = sys.argv[2]
cohort_prefix = sys.argv[3]
bucket_id = sys.argv[4]

storage_client = storage.Client(project='talkowski-sv-gnomad')
bucket = storage_client.get_bucket(bucket_id.replace('gs://', ''))

sample_table = pd.read_csv(sample_tsv_uri, sep='\t')
ped = pd.read_csv(ped_uri, sep='\t')
ped.columns = ['FamID', 'IndividualID', 'FatherID', 'MotherID', 'Gender', 'Affected']

sample_id_dict = sample_table.set_index('entity:sample_id').to_dict()['collaborator_id']
sample_id_dict['0'] = '0'

for sample_id in ['IndividualID', 'FatherID', 'MotherID']:
    ped[sample_id] = ped[sample_id].astype(str).map(sample_id_dict)

ped.to_csv(f"{cohort_prefix}.ped", sep='\t', index=False)
bucket.blob(f"resources/pedigrees/{cohort_prefix}.ped").upload_from_string(ped.to_csv(sep='\t', index=False), 'text/csv')

