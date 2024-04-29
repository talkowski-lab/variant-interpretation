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

workflow step3 {
    input {
        Array[String] filtered_mt
        File ped_sex_qc
        File loeuf_file
        Boolean hail_autoscale
        String bucket_id
        String cohort_prefix
        String hail_denovo_filtering_script
        String hail_docker
        String sv_base_mini_docker
        Float min_child_ab=0.25
        Float min_dp_ratio=0.1
        Float min_gq=25
    }

    scatter (mt_uri in filtered_mt) {
        call helpers.getHailMTSize as getStep2MTSize {
            input:
                mt_uri=mt_uri,
                hail_docker=hail_docker
        }
    
        call hailDenovoFilteringRemote {
            input:
                filtered_mt=mt_uri,
                input_size=getStep2MTSize.mt_size,
                ped_sex_qc=ped_sex_qc,
                bucket_id=bucket_id,
                cohort_prefix=cohort_prefix,
                loeuf_file=loeuf_file,
                hail_denovo_filtering_script=hail_denovo_filtering_script,
                hail_docker=hail_docker
        }
    }

    output {
        # step 3 output
        Array[File] de_novo_results_sharded = hailDenovoFilteringRemote.de_novo_results
        Array[File] de_novo_vep_sharded = hailDenovoFilteringRemote.de_novo_vep
        Array[String] de_novo_ht = hailDenovoFilteringRemote.de_novo_ht
        Array[String] tdt_mt = hailDenovoFilteringRemote.tdt_mt
        Array[String] tdt_parent_aware_mt = hailDenovoFilteringRemote.tdt_parent_aware_mt
    }
}


task hailDenovoFilteringRemote {
    input {
        File ped_sex_qc
        Float input_size
        String filtered_mt
        String bucket_id
        String cohort_prefix
        String loeuf_file
        String hail_denovo_filtering_script
        String hail_docker
        Float min_child_ab
        Float min_dp_ratio
        Float min_gq
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
        curl ~{hail_denovo_filtering_script} > hail_denovo_filtering_script.py
        python3 hail_denovo_filtering_script.py ~{filtered_mt} ~{cohort_prefix} ~{ped_sex_qc} ~{loeuf_file} \
        ~{cpu_cores} ~{memory} ~{bucket_id} ~{min_child_ab} ~{min_dp_ratio} ~{min_gq} > stdout
    }

    String prefix = basename(filtered_mt, "_wes_denovo_basic_filtering.mt")
    output {
        File de_novo_results = "~{prefix}_wes_final_denovo.txt"
        File de_novo_vep = "~{prefix}_wes_final_denovo_vep.txt"
        String de_novo_ht = read_lines('mt_uri.txt')[0]
        String tdt_mt = read_lines('mt_uri.txt')[1]
        String tdt_parent_aware_mt = read_lines('mt_uri.txt')[2]
    }
}
