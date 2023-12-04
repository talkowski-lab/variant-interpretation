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
    python3 ~{convert_ped_sampleIDs_python_script} ~{sample_tsv_uri} ~{ped_uri} ~{cohort_prefix} ~{bucket_id}
    >>>

    output {
        File ped_uri = "~{bucket_id}/resources/pedigrees/~{cohort_prefix}.ped"
    }
}