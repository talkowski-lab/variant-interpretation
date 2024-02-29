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
        String bucket_id
        String hail_docker
    }

    call getBucketSizes {
        input:
        bucket_id=bucket_id,
        hail_docker=hail_docker
    }

    output {
        File file_sizes = getBucketSizes.file_sizes
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
        
        bucket_id = sys.argv[1]

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