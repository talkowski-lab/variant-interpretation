from typing import Union
from tenacity import retry, after_log, before_sleep_log, retry_if_exception_type, stop_after_attempt, wait_exponential
from firecloud import fiss
from firecloud.errors import FireCloudServerError
import firecloud.api as fapi
import pandas as pd 
import numpy as np
import re
import sys
import time
import math
import ast
import os
import datetime

import logging
from logging import INFO, DEBUG
logger = logging.getLogger()
logger.setLevel(INFO)

from google.cloud import storage
import pandas as pd 
import numpy as np
import datetime
import sys
from terra_data_table_util import get_terra_table_to_df

BILLING_PROJECT_ID = sys.argv[1]
WORKSPACE = sys.argv[2]
BUCKET = sys.argv[3]

bucket_id = BUCKET.split('gs://')[1]
storage_client = storage.Client()
bucket = storage_client.bucket(bucket_id)

cohort_df = get_terra_table_to_df(project=BILLING_PROJECT_ID, workspace=WORKSPACE, table_name='cohort')
cohort_df.index = cohort_df['entity:cohort_id']

suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
def humansize(nbytes):
    i = 0
    while nbytes >= 1024 and i < len(suffixes)-1:
        nbytes /= 1024.
        i += 1
    f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
    return '%s %s' % (f, suffixes[i])

submissions = pd.DataFrame(fapi.list_submissions(BILLING_PROJECT_ID, WORKSPACE).json()).sort_values('submissionDate',ascending=False).set_index('submissionId',drop=False)
submissions = pd.concat([submissions, submissions.submissionEntity.apply(pd.Series)], axis=1)
submissions['submissionDate'] = submissions.submissionDate.apply(lambda date_str: datetime.datetime.fromisoformat(date_str[:-1]))

# Get output names for each method/task
submissions['methodConfigurationNameClean'] = submissions.methodConfigurationName.str.split('_').str[0]
method_outputs = {}
for x in submissions.methodConfigurationNameClean.unique():
    try:
        method_outputs[x] = [x.split('.')[1] for x in fapi.get_workspace_config(config=x, 
                          namespace=BILLING_PROJECT_ID, 
                          cnamespace=BILLING_PROJECT_ID, workspace=WORKSPACE).json()['outputs'].values()]
    except:
        method_outputs[x] = np.nan

only_successful = submissions[submissions.workflowStatuses.apply(lambda dct: 'Succeeded' in list(dct.keys()))]
only_failed = submissions[~submissions.index.isin(only_successful.index)]

# Successful submissions first
tot_subs = only_successful.shape[0]
i = 0
successful_submissions = pd.DataFrame()

# Get workflow info
print("------ Getting workflow info for each submission ------")
for sub_id, submission in only_successful.iterrows():
    i += 1
    if (i%10==0): 
        print(f"submission {i}/{tot_subs}...")
    metadata = pd.DataFrame(fapi.get_submission(submission_id=sub_id, namespace=BILLING_PROJECT_ID, workspace=WORKSPACE).json()['workflows'])
    try:
        submission['workflowId'] = metadata['workflowId'].tolist()
        submission['workflowStatus'] = metadata['status'].tolist()
        submission['workflowEntity'] = metadata['workflowEntity'].tolist()
    except:
        pass
    successful_submissions = pd.concat([successful_submissions, pd.DataFrame(submission).T])

# Get cohort sets into separate rows by cohort
successful_submissions_exploded = successful_submissions.explode(['workflowId','workflowStatus','workflowEntity'])
# update entityName and entityType
successful_submissions_exploded = successful_submissions_exploded.workflowEntity.apply(pd.Series).combine_first(successful_submissions_exploded)

successful_submissions_exploded['method_outputs'] = successful_submissions_exploded.methodConfigurationNameClean.map(method_outputs)

# Get actual outputs for each submission and workflow
print("------ Getting workflow outputs for each submission ------")
successful_submission_outputs = pd.DataFrame()
tot_subs = successful_submissions_exploded.shape[0]
i = 0
for sub_id, submission in successful_submissions_exploded.iterrows():
    i += 1
    if (i%10==0): 
        print(f"submission {i}/{tot_subs}...")
    outputs = {k.split('.')[1]: v for k, v in fapi.get_workflow_metadata(namespace=BILLING_PROJECT_ID, workspace=WORKSPACE, submission_id=sub_id, workflow_id=submission.workflowId).json()['outputs'].items()}
    submission['method_outputs'] = outputs.keys()
    submission['workflow_outputs'] = outputs.values()
    successful_submission_outputs = pd.concat([successful_submission_outputs, pd.DataFrame(submission).T])    

# Get multiple outputs for a method into separate rows
successful_submission_outputs = successful_submission_outputs.explode(['method_outputs', 'workflow_outputs'])

# Check if in Terra data table
def in_data_table(submission):
    cohort = submission.entityName
    output_name = submission.method_outputs
    try:
        workflow_output = submission.workflow_outputs
        if '[' in workflow_output:
            workflow_output = ast.literal_eval(workflow_output)
        data_table_val = cohort_df.loc[cohort, output_name]
        if '[' in data_table_val:
            data_table_val = ast.literal_eval(data_table_val)
        return data_table_val==workflow_output
    except Exception as e:
        return np.nan

successful_submission_outputs['in_data_table'] = successful_submission_outputs.apply(in_data_table, axis=1)

# Sanity check that we're keeping one of each output for each cohort
entity_outputs_to_remove = successful_submission_outputs[successful_submission_outputs.in_data_table==False][['entityName','method_outputs']].value_counts().reset_index()
entity_outputs_to_keep = successful_submission_outputs[successful_submission_outputs.in_data_table==True][['entityName','method_outputs']].value_counts().reset_index()

entity_outputs_to_remove['entityName:method_outputs'] = entity_outputs_to_remove[['entityName','method_outputs']].agg(':'.join, axis=1)
entity_outputs_to_keep['entityName:method_outputs'] = entity_outputs_to_keep[['entityName','method_outputs']].agg(':'.join, axis=1)

print("--- Check that all entities and method outputs being removed have an updated version ---")
print(entity_outputs_to_remove['entityName:method_outputs'].isin(entity_outputs_to_keep['entityName:method_outputs']).value_counts())

n_filled_in_data_table = (~cohort_df.isna()).sum()[np.intersect1d(cohort_df.columns, successful_submission_outputs.method_outputs.dropna().unique())].sum()
print(f"Terra data table has {n_filled_in_data_table} total entries filled.")

# Get sizes
def get_blob_size_directory(uri):
    tot_size = 0
    blobs_list = storage_client.list_blobs(bucket, prefix=uri.split(BUCKET)[1][1:])
    for blob in blobs_list:
        tot_size += blob.size
    return tot_size

def get_blob_size(workflow_outputs):
    try:
        if type(workflow_outputs)!=list:
            workflow_outputs = [workflow_outputs]
        tot_size = 0
        for uri in workflow_outputs:
            if '.mt' in uri or '.ht' in uri:
                tot_size += get_blob_size_directory(uri)
            else:
                blob = bucket.get_blob(uri.split(BUCKET)[1][1:])
                tot_size += blob.size
        return tot_size
    except Exception as e:
        return np.nan
    
print("--- Getting file sizes for successful submissions' outputs ---")
successful_submission_outputs['sizes'] = successful_submission_outputs.workflow_outputs.apply(get_blob_size)

# Now get sizes for failed submission directories
print("--- Getting file sizes for failed submissions' outputs ---")
only_failed['uri'] = BUCKET + '/submissions/' + only_failed.submissionId + '/'
only_failed['sizes'] = only_failed.uri.apply(get_blob_size_directory)

successful_submission_outputs.to_csv(f"{str(datetime.datetime.now().strftime('%Y-%m-%d'))}_successful_submissions.tsv", sep='\t', index=False)
only_failed.to_csv(f"{str(datetime.datetime.now().strftime('%Y-%m-%d'))}_failed_submissions.tsv", sep='\t', index=False)

def split_df_by_column_sum(df, column_name, target_sum):
    current_sum = 0
    current_group = 0
    new_df = pd.DataFrame()
    for sub_id, submission in df.iterrows():
        current_sum += submission[column_name]
        submission['group'] = current_group
        if current_sum >= target_sum:
            current_group += 1
            current_sum = 0
        new_df = pd.concat([new_df, pd.DataFrame(submission).T])
    return new_df

# Chunk failed submissions
chunk_size = 1_000_000_000  # 1 TB
failed_submissions_to_delete = split_df_by_column_sum(only_failed, 'sizes', chunk_size)
failed_submissions_to_delete['group'] = failed_submissions_to_delete.group.astype(int)
for group in failed_submissions_to_delete.group.unique():
    group_df = failed_submissions_to_delete[failed_submissions_to_delete.group==group]
    to_delete_filename = f"failed_terra_submissions_to_delete_chunk_{group}_{str(datetime.datetime.now().strftime('%Y-%m-%d'))}.txt"
    group_df.uri.to_csv(to_delete_filename, index=False, header=None)  # these are directories!

# Chunk successful submissions
successful_submissions_to_delete = successful_submission_outputs[successful_submission_outputs.in_data_table==False]
successful_submissions_to_delete = split_df_by_column_sum(successful_submissions_to_delete, 'sizes', chunk_size)
successful_submissions_to_delete['group'] = successful_submissions_to_delete.group.astype(int)
for group in successful_submissions_to_delete.group.unique():
    group_df = successful_submissions_to_delete[successful_submissions_to_delete.group==group]
    to_delete_filename = f"terra_submissions_to_delete_chunk_{group}_{str(datetime.datetime.now().strftime('%Y-%m-%d'))}.txt"
    group_outputs = [file if type(lst)==list else lst for lst in group_df.workflow_outputs.tolist() for file in lst]
    pd.Series(group_outputs).to_csv(to_delete_filename, index=False, header=None)

space_cleared = humansize(failed_submissions_to_delete.sizes.sum() + successful_submissions_to_delete.sizes.sum())
with open("space_cleared.txt", "w") as text_file:
    print(space_cleared, file=text_file)
