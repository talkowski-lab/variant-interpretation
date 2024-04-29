version 1.0

import "wes-denovo-helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow step2 {
    input {
        File ped_sex_qc
        File lcr_uri
        Array[String] annot_mt
        String cohort_prefix
        String hail_basic_filtering_script
        String hail_docker
        String bucket_id
        Float call_rate_threshold
        RuntimeAttr? runtime_attr_override
    }

    scatter (mt_uri in annot_mt) {
        call helpers.getHailMTSize as getHailMTSize {
            input:
                mt_uri=mt_uri,
                hail_docker=hail_docker
        }
        call hailBasicFilteringRemote {
            input:
                lcr_uri=lcr_uri,
                annot_mt=mt_uri,
                input_size=select_first([getHailMTSize.mt_size]),
                ped_sex_qc=ped_sex_qc,
                bucket_id=bucket_id,
                cohort_prefix=cohort_prefix,
                hail_basic_filtering_script=hail_basic_filtering_script,
                hail_docker=hail_docker,
                call_rate_threshold=call_rate_threshold,
                runtime_attr_override=runtime_attr_override
        }
    }

    output {
        Array[String] filtered_mt = hailBasicFilteringRemote.filtered_mt
        Array[File] post_filter_sample_qc_info = hailBasicFilteringRemote.post_filter_sample_qc_info
    }
}

task hailBasicFilteringRemote {
    input {
        File ped_sex_qc
        File lcr_uri
        Float input_size
        Float call_rate_threshold
        String annot_mt
        String bucket_id
        String cohort_prefix
        String hail_basic_filtering_script
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
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

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command {
        curl ~{hail_basic_filtering_script} > hail_basic_filtering_script.py
        python3 hail_basic_filtering_script.py ~{annot_mt} ~{cohort_prefix} ~{ped_sex_qc} \
        ~{cpu_cores} ~{memory} ~{bucket_id} ~{lcr_uri} ~{call_rate_threshold}
    }

    String prefix = basename(annot_mt, "_wes_denovo_annot.mt")
    output {
        String filtered_mt = read_lines("mt_uri.txt")[0]
        File post_filter_sample_qc_info = "~{prefix}_wes_final_annot_post_filter_qc_info.txt"
    }
}
