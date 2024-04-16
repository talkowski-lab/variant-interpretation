version 1.0 

import "prioritizeCSQ.wdl" as prioritizeCSQ

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow step6 {
    input {
        File vcf_metrics_tsv_annot
        Float AF_threshold=0.005
        Int AC_threshold=2
        Float csq_af_threshold=0.01
        String filter_final_tsv_script
        String prioritize_csq_script
        String hail_docker
        String sample_column
        RuntimeAttr? runtime_attr_prioritize
        RuntimeAttr? runtime_attr_filter_final
    }

    call prioritizeCSQ.annotateMostSevereCSQ as prioritizeCSQ {
        input:
        vcf_metrics_tsv=vcf_metrics_tsv_annot,
        prioritize_csq_script=prioritize_csq_script,
        hail_docker=hail_docker,
        sample_column=sample_column,
        runtime_attr_override=runtime_attr_prioritize
    }

    call filterFinalTSV {
        input:
            vcf_metrics_tsv_annot=prioritizeCSQ.vcf_metrics_tsv_prior_csq,
            AF_threshold=AF_threshold,
            AC_threshold=AC_threshold,
            csq_af_threshold=csq_af_threshold,
            filter_final_tsv_script=filter_final_tsv_script,
            hail_docker=hail_docker,
            runtime_attr_filter_final=runtime_attr_filter_final
    }

    output {
        File vcf_metrics_tsv_final = filterFinalTSV.vcf_metrics_tsv_final
    }
}

task filterFinalTSV {
    input {
        File vcf_metrics_tsv_annot
        Float AF_threshold
        Int AC_threshold
        Float csq_af_threshold
        String filter_final_tsv_script
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_metrics_tsv_annot, "GB")
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
        curl ~{filter_final_tsv_script} > filter_tsv.py
        python3 filter_tsv.py ~{vcf_metrics_tsv_annot} ~{AC_threshold} ~{AF_threshold} ~{csq_af_threshold} 
    }

    output {
        File vcf_metrics_tsv_final = basename(vcf_metrics_tsv_annot, '.tsv.gz') + '_filtered.tsv.gz'
    }
}