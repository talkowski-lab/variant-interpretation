version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow getTerraBucketSizes {
    input {
        String BILLING_PROJECT_ID
        String WORKSPACE
        String bucket_id
        String hail_docker
    }

    call getBucketSizes {
        input:
        bucket_id=bucket_id,
        hail_docker=hail_docker
    }

    call getSubmissionInfo {
        input:
        file_sizes=getBucketSizes.file_sizes,
        BILLING_PROJECT_ID=BILLING_PROJECT_ID,
        WORKSPACE=WORKSPACE,
        bucket_id=bucket_id,
        hail_docker=hail_docker
    }

    call getSubmissionsToDelete {
        input:
        file_sizes=getBucketSizes.file_sizes,
        submission_info=getSubmissionInfo.submission_info,
        BILLING_PROJECT_ID=BILLING_PROJECT_ID,
        WORKSPACE=WORKSPACE,
        bucket_id=bucket_id,
        hail_docker=hail_docker
    }

    call deleteSubmissions {
        input:
        submissions_to_delete=getSubmissionsToDelete.submissions_to_delete,
        hail_docker=hail_docker
    }

    output {
        File file_sizes = getBucketSizes.file_sizes
        File submission_info = getSubmissionInfo.submission_info
        File submissions_to_delete = getSubmissionsToDelete.submissions_to_delete
        String space_cleared = getSubmissionsToDelete.space_cleared
    }
}

task getBucketSizes {
    input {
        String bucket_id
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float base_disk_gb = 10.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
    
    command <<<
        cat <<EOF > get_bucket_sizes.py 
        from google.cloud import storage
        import pandas as pd 
        import numpy as np
        import datetime
        import sys
        
        bucket_id = sys.argv[1].strip('gs://')

        storage_client = storage.Client()
        blobs_list = storage_client.list_blobs(bucket_or_name=bucket_id)

        blob_ids = []
        blob_sizes = []
        for blob in blobs_list:
            if blob.size > 1_000_000:  # size is in bytes
                blob_ids.append(blob.id)
                blob_sizes.append(blob.size)

        suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
        def humansize(nbytes):
            i = 0
            while nbytes >= 1024 and i < len(suffixes)-1:
                nbytes /= 1024.
                i += 1
            f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
            return '%s %s' % (f, suffixes[i])
        large_files = pd.DataFrame({'uri': blob_ids, 'sizes': blob_sizes,
                           'human_sizes': [humansize(s) for s in blob_sizes]}).sort_values(by='sizes', ascending=False).reset_index(drop=True)

        large_files.to_csv(f"terra_storage_file_sizes_{str(datetime.datetime.now().strftime('%Y-%m-%d'))}.tsv", sep='\t', index=False)
        EOF

        python3 get_bucket_sizes.py ~{bucket_id} 
    >>>

    output {
        File file_sizes = glob('terra_storage_file_sizes*')[0]
    }
}

task getSubmissionInfo {
    input {
        File file_sizes
        String BILLING_PROJECT_ID
        String WORKSPACE
        String bucket_id
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(file_sizes, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    
    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
    cat <<EOF > get_submission_info.py
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

    BILLING_PROJECT_ID = sys.argv[1]
    WORKSPACE = sys.argv[2]
    BUCKET = sys.argv[3]
    file_sizes = sys.argv[4]

    large_files = pd.read_csv(file_sizes, sep='\t')

    suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
    def humansize(nbytes):
        i = 0
        while nbytes >= 1024 and i < len(suffixes)-1:
            nbytes /= 1024.
            i += 1
        f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
        return '%s %s' % (f, suffixes[i])

    submissions = pd.DataFrame(fapi.list_submissions(BILLING_PROJECT_ID, WORKSPACE).json()).sort_values('submissionDate',ascending=False)
    submissions = pd.concat([submissions, submissions.submissionEntity.apply(pd.Series)], axis=1)

    submission_files = large_files[large_files.uri.str.split('/').str[1]=='submissions'].copy()
    submission_files.loc[:,'submission_id'] = submission_files.uri.str.split('/').str[2]

    submission_id_sizes = pd.DataFrame(submission_files.groupby('submission_id').sizes.sum()).reset_index().sort_values(by='sizes', ascending=False)
    submission_id_sizes['human_sizes'] = submission_id_sizes.sizes.apply(humansize)
    submission_id_sizes['status'] = submission_id_sizes.submission_id.map(submissions.set_index('submissionId').status.to_dict())
    submission_id_sizes['submissionDate'] = submission_id_sizes.submission_id.map(submissions.set_index('submissionId').submissionDate.to_dict())

    def get_workflow_status(sub_id):
        submission = fapi.get_submission(BILLING_PROJECT_ID, WORKSPACE, sub_id).json()

        workflow_id = submission['workflows'][0]['workflowId']
        try:
            workflow = fapi.get_workflow_metadata(BILLING_PROJECT_ID, WORKSPACE, sub_id, workflow_id).json()
            return workflow['status']
        except:
            return np.nan

    submission_id_sizes['workflow_status'] = submission_id_sizes.submission_id.apply(get_workflow_status)
    
    submission_id_sizes.to_csv(f"terra_submissions_{str(datetime.datetime.now().strftime('%Y-%m-%d'))}.tsv", sep='\t', index=False)
    EOF

    python3 get_submission_info.py ~{BILLING_PROJECT_ID} ~{WORKSPACE} ~[bucket_id] ~{file_sizes}
    >>>

    output {
        File submission_info = glob("terra_submissions_*")[0]
    }
}

task getSubmissionsToDelete {
    input {
        File file_sizes
        File submission_info
        String BILLING_PROJECT_ID
        String WORKSPACE
        String bucket_id
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size([file_sizes, submission_info], "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    
    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
    cat <<EOF > get_submissions_to_delete.py
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

    BILLING_PROJECT_ID = sys.argv[1]
    WORKSPACE = sys.argv[2]
    BUCKET = sys.argv[3]
    file_sizes = sys.argv[4]
    submission_info = sys.argv[5]

    suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
    def humansize(nbytes):
        i = 0
        while nbytes >= 1024 and i < len(suffixes)-1:
            nbytes /= 1024.
            i += 1
        f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
        return '%s %s' % (f, suffixes[i])

    large_files = pd.read_csv(file_sizes, sep='\t')

    large_files.loc[:,'submission_id'] = large_files.uri.str.split('/').str[2]

    submission_info = pd.read_csv(submission_info, sep='\t').set_index('submission_id', drop=False)

    submission_info['submissionDate'] = submission_info.submissionDate.apply(lambda date_str: datetime.datetime.fromisoformat(date_str[:-1]))
    submission_info['submission_path'] = BUCKET + '/submissions/' + submission_info.submission_id + '/'

    today = datetime.datetime.today()
    week_ago = today - datetime.timedelta(days=7)

    ## All old failed submissions (low-hanging fruit)
    # all failed/aborted submissions from over a week ago -- delete
    submissions_to_delete = submission_info[(submission_info.submissionDate < week_ago)
                                            & (submission_info.workflow_status!='Succeeded')].submission_id
    failed_new_submissions = submission_info[(submission_info.submissionDate >= week_ago)
                                            & (submission_info.workflow_status!='Succeeded')].submission_id

    ## Old successful submissions that have been overrided
    old_successful_submissions = submission_info[(submission_info.workflow_status=='Succeeded')].submission_id

    submissions = pd.DataFrame(fapi.list_submissions(BILLING_PROJECT_ID, WORKSPACE).json()).sort_values('submissionDate',ascending=False).set_index('submissionId',drop=False)
    submissions = pd.concat([submissions, submissions.submissionEntity.apply(pd.Series)], axis=1)
    submissions['submissionDate'] = submissions.submissionDate.apply(lambda date_str: datetime.datetime.fromisoformat(date_str[:-1]))
    submissions['workflow_status'] = submissions.submissionId.map(submission_info.workflow_status.to_dict())

    def get_all_submissions(submission_id):
        workflow_name = submissions[submissions.submissionId==submission_id].methodConfigurationName.tolist()[0].strip('_')[0]
        all_workflow_submissions = submissions[submissions.methodConfigurationName.str.contains(workflow_name)]
        return all_workflow_submissions

    def has_updated_submission(submission_id):
        workflow_name = submissions[submissions.submissionId==submission_id].methodConfigurationName.tolist()[0].split('_')[0]
        all_workflow_submissions = submissions[(submissions.methodConfigurationName.str.contains(workflow_name))]
        try:
            cohort = submissions[submissions.submissionId==submission_id].submissionEntity[0]['entityName']
            all_workflow_submissions = all_workflow_submissions[all_workflow_submissions.submissionEntity.apply(lambda entity_dict:
                                                                                                            entity_dict['entityName']==cohort)]
        except:
            pass
        successful_sorted_runs = all_workflow_submissions[all_workflow_submissions.workflow_status=='Succeeded'].sort_values('submissionDate', ascending=False)
        return successful_sorted_runs.iloc[0,:].submissionId==submission_id

    old_submission_df = pd.DataFrame(submission_info[submission_info.workflow_status=='Succeeded'].submission_id)
    old_submission_df['is_most_recent'] = old_submission_df.submission_id.apply(has_updated_submission)

    old_successful_submissions_to_delete = old_submission_df[~old_submission_df.is_most_recent].submission_id
    all_submissions_to_delete = pd.concat([submissions_to_delete, failed_new_submissions, old_successful_submissions_to_delete])
    paths_to_delete = submission_info.loc[all_submissions_to_delete].submission_path

    to_delete_filename = f"terra_submissions_to_delete_{str(datetime.datetime.now().strftime('%Y-%m-%d'))}.txt"
    paths_to_delete.to_csv(to_delete_filename, index=False, header=None)
    space_cleared = humansize(large_files[large_files.submission_id.isin(all_submissions_to_delete)].sizes.sum())
    with open("space_cleared.txt", "w") as text_file:
        print(space_cleared, file=text_file)
    EOF

    python3 get_submissions_to_delete.py ~{BILLING_PROJECT_ID} ~{WORKSPACE} ~{bucket_id} ~{file_sizes} ~{submission_info} 
    >>>

    output {
        File submissions_to_delete = glob("terra_submissions_to_delete_*")[0]
        String space_cleared = read_lines("space_cleared.txt")[0]
    }
}

task deleteSubmissions {
    input {
        File submissions_to_delete
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(submissions_to_delete, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    
    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command {
        set -eou pipefail
        gsutil -m rm -r $(cat ~{submissions_to_delete})
    }
}