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

workflow WESprioritizeCSQ {
    input {
        File de_novo_results
        File de_novo_vep
        String cohort_prefix
        String prioritize_csq_script
        String hail_docker
        String sample_column
    }

    call mergeVEPIntoResults {
        input:
        de_novo_results=de_novo_results,
        de_novo_vep=de_novo_vep,
        cohort_prefix=cohort_prefix,
        hail_docker=hail_docker
    }

    call prioritizeCSQ.annotateMostSevereCSQ as prioritizeCSQ {
        input:
        vcf_metrics_tsv=mergeVEPIntoResults.de_novo_merged,
        prioritize_csq_script=prioritize_csq_script,
        hail_docker=hail_docker,
        sample_column=sample_column
    }

    output {
        File de_novo_merged = prioritizeCSQ.vcf_metrics_tsv_prior_csq
    }
}

task mergeVEPIntoResults {
    input {
        File de_novo_results
        File de_novo_vep
        String cohort_prefix
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(de_novo_results, "GB") + size(de_novo_vep, "GB")
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
    cat <<EOF > merge_vep.py
    import pandas as pd
    import sys

    de_novo_results = sys.argv[1]
    de_novo_vep = sys.argv[2]
    cohort_prefix = sys.argv[3]

    dn_merged_df = pd.read_csv(de_novo_results, sep='\t')
    dn_vep_df = pd.read_csv(de_novo_vep, sep='\t')
    dn_merged_df.merge(dn_vep_df).to_csv(f"{cohort_prefix}_wes_final_denovo_vep_merged.tsv", sep='\t', index=False)

    EOF

    python3 merge_vep.py ~{de_novo_results} ~{de_novo_vep} ~{cohort_prefix}
    >>>

    output {
        File de_novo_merged = "~{cohort_prefix}_wes_final_denovo_vep_merged.tsv"
    }
}
