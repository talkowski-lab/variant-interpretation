version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow getTerraSubmissionSizes {
    input {
        String BILLING_PROJECT_ID
        String WORKSPACE
        String bucket_id
        String hail_docker
        String terra_data_table_util_script
        String get_submission_sizes_script
        Boolean delete_submissions=false
    }

    call getSubmissionSizes {
        input:
        BILLING_PROJECT_ID=BILLING_PROJECT_ID,
        WORKSPACE=WORKSPACE,
        bucket_id=bucket_id,
        hail_docker=hail_docker,
        terra_data_table_util_script=terra_data_table_util_script,
        get_submission_sizes_script=get_submission_sizes_script
    }

    # if (delete_submissions) {
    #     scatter (submissions_shard in getSubmissionsToDelete.submissions_to_delete) {
    #         call deleteSubmissions {
    #             input:
    #             submissions_to_delete=submissions_shard,
    #             hail_docker=hail_docker
    #         }
    #     }
    # }
    output {
        File successful_submissions = getSubmissionSizes.successful_submissions
        File failed_submissions = getSubmissionSizes.failed_submissions
        # Array[File] submissions_to_delete = getSubmissionsToDelete.submissions_to_delete
        # String space_cleared = getSubmissionsToDelete.space_cleared
    }
}

task getSubmissionSizes {
    input {
        String BILLING_PROJECT_ID
        String WORKSPACE
        String bucket_id
        String hail_docker
        String terra_data_table_util_script
        String get_submission_sizes_script
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
        curl ~{terra_data_table_util_script} > terra_data_table_util.py
        curl ~{get_submission_sizes_script} > get_submission_sizes.py
        python3 get_submission_sizes.py ~{BILLING_PROJECT_ID} ~{WORKSPACE} ~{bucket_id} > stdout
    >>>

    output {
        File successful_submissions = glob('_successful_submissions.tsv')[0]
        File failed_submissions = glob('_failed_submissions.tsv')[0]
    }
}