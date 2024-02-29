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

    output {
        File file_sizes = getBucketSizes.file_sizes
        File submission_info = getSubmissionInfo.submission_info
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
        workflow = fapi.get_workflow_metadata(BILLING_PROJECT_ID, WORKSPACE, sub_id, workflow_id).json()

        return workflow['status']

    submission_id_sizes['workflow_status'] = submission_id_sizes.submission_id.apply(get_workflow_status)
    
    submission_id_sizes.to_csv(f"terra_submissions_{str(datetime.datetime.now().strftime('%Y-%m-%d'))}.tsv", sep='\t', index=False)
    EOF

    python3 get_submission_info.py ~{BILLING_PROJECT_ID} ~{WORKSPACE} ~[bucket_id] ~{file_sizes}
    >>>

    output {
        File submission_info = glob("terra_submissions_*")[0]
    }
}