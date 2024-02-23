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
        File ped_uri
        File filtered_mt
        File loeuf_file
        String cohort_prefix
        String hail_denovo_filtering_script
        String hail_docker
        String bucket_id
        RuntimeAttr? runtime_attr_override
    }

    if (bucket_id=='false') {
        call hailDenovoFiltering {
            input:
                filtered_mt=filtered_mt,
                ped_uri=ped_uri,
                cohort_prefix=cohort_prefix,
                loeuf_file=loeuf_file,
                hail_denovo_filtering_script=hail_denovo_filtering_script,
                hail_docker=hail_docker,
                runtime_attr_override=runtime_attr_override
        }
    }

    if (bucket_id!='false') {
        call helpers.getHailMTSize as getHailMTSize {
            input:
                mt_uri=filtered_mt,
                hail_docker=hail_docker
        }
        call hailDenovoFilteringRemote {
            input:
                input_size=select_first([getHailMTSize.mt_size]),
                filtered_mt=filtered_mt,
                ped_uri=ped_uri,
                bucket_id=bucket_id,
                cohort_prefix=cohort_prefix,
                loeuf_file=loeuf_file,
                hail_denovo_filtering_script=hail_denovo_filtering_script,
                hail_docker=hail_docker,
                runtime_attr_override=runtime_attr_override
        }
    }

    output {
        File de_novo_results = select_first([hailDenovoFiltering.de_novo_results, hailDenovoFilteringRemote.de_novo_results])
        File de_novo_vep = select_first([hailDenovoFiltering.de_novo_vep, hailDenovoFilteringRemote.de_novo_vep])
        String de_novo_ht = select_first([hailDenovoFiltering.de_novo_ht, hailDenovoFilteringRemote.de_novo_ht])
        String tdt_mt = select_first([hailDenovoFiltering.tdt_mt, hailDenovoFilteringRemote.tdt_mt])
        String tdt_parent_aware_mt = select_first([hailDenovoFiltering.tdt_parent_aware_mt, hailDenovoFilteringRemote.tdt_parent_aware_mt])
    }
}

task hailDenovoFiltering {
    input {
        File ped_uri
        File filtered_mt
        String cohort_prefix
        String loeuf_file
        String hail_denovo_filtering_script
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(filtered_mt, "GB")
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

    Boolean bucket_id = false

    command {
        tar -zxvf ~{filtered_mt}

        curl ~{hail_denovo_filtering_script} > hail_denovo_filtering_script.py
        python3 hail_denovo_filtering_script.py ~{basename(filtered_mt, '.gz')} ~{cohort_prefix} ~{ped_uri} ~{loeuf_file} \
        ~{cpu_cores} ~{memory} ~{bucket_id}

        tar -zcvf ~{cohort_prefix}_wes_final_denovo.ht.gz ~{cohort_prefix}_wes_final_denovo.ht
        tar -zcvf ~{cohort_prefix}_trio_tdt.mt.gz ~{cohort_prefix}_trio_tdt.mt
        tar -zcvf ~{cohort_prefix}_parent_aware_trio_tdt.mt.gz ~{cohort_prefix}_parent_aware_trio_tdt.mt
    }

    output {
        File de_novo_results = "~{cohort_prefix}_wes_final_denovo.txt"
        File de_novo_vep = "~{cohort_prefix}_wes_final_denovo_vep.txt"
        File de_novo_ht = "~{cohort_prefix}_wes_final_denovo.ht.gz"
        File tdt_mt = "~{cohort_prefix}_trio_tdt.mt.gz"
        File tdt_parent_aware_mt = "~{cohort_prefix}_parent_aware_trio_tdt.mt.gz"
    }
}

task hailDenovoFilteringRemote {
    input {
        File ped_uri
        Float input_size
        String filtered_mt
        String bucket_id
        String cohort_prefix
        String loeuf_file
        String hail_denovo_filtering_script
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
        curl ~{hail_denovo_filtering_script} > hail_denovo_filtering_script.py
        python3 hail_denovo_filtering_script.py ~{filtered_mt} ~{cohort_prefix} ~{ped_uri} ~{loeuf_file} \
        ~{cpu_cores} ~{memory} ~{bucket_id}
    }

    output {
        File de_novo_results = "~{cohort_prefix}_wes_final_denovo.txt"
        File de_novo_vep = "~{cohort_prefix}_wes_final_denovo_vep.txt"
        String de_novo_ht = read_lines('mt_uri.txt')[0]
        String tdt_mt = read_lines('mt_uri.txt')[1]
        String tdt_parent_aware_mt = read_lines('mt_uri.txt')[2]
    }
}
