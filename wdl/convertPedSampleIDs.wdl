version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow convertPedSampleIDs {
    input {
        File sample_tsv_uri
        File raw_ped
        String convert_ped_sampleIDs_python_script
        String cohort_prefix
        String bucket_id
        String hail_docker
    }

    call runConvertPed {
        input:
            sample_tsv_uri=sample_tsv_uri,
            raw_ped=raw_ped,
            convert_ped_sampleIDs_python_script=convert_ped_sampleIDs_python_script,
            cohort_prefix=cohort_prefix,
            bucket_id=bucket_id,
            hail_docker=hail_docker
    }

    output {
        File ped_uri = runConvertPed.ped_uri
    }
}

task runConvertPed {
    input {
        File sample_tsv_uri
        File raw_ped
        String convert_ped_sampleIDs_python_script
        String cohort_prefix
        String bucket_id
        String hail_docker
    }

    runtime {
        docker: hail_docker
    }

    command <<<
    curl ~{convert_ped_sampleIDs_python_script} > convert.py
    python3 convert.py ~{sample_tsv_uri} ~{raw_ped} ~{cohort_prefix} ~{bucket_id}
    >>>

    output {
        File ped_uri = cohort_prefix + ".ped"
    }
}