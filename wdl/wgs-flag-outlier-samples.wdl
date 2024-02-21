version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow flagOutlierSamples {
    input {
        File vcf_metrics_tsv
        String cohort_prefix
        String hail_docker
    }

    call getOutliers {
        input:
        vcf_metrics_tsv=vcf_metrics_tsv,
        cohort_prefix=cohort_prefix,
        hail_docker=hail_docker
    }

    output {
        Array[String] outlier_samples = getOutliers.outlier_samples
    }
}

task getOutliers {
    input {
        File vcf_metrics_tsv
        String cohort_prefix
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_metrics_tsv, "GB")
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

    command <<<
        cat <<EOF > get_outliers.py 
        import os
        import sys
        import pandas as pd
        import numpy as np

        cohort_prefix = sys.argv[1]
        vcf_metrics_tsv = sys.argv[2]

        df = pd.read_csv(vcf_metrics_tsv, sep='\t')
        snv_counts = df[df.TYPE=='SNV'].SAMPLE.value_counts()
        indel_counts = df[df.TYPE=='Indel'].SAMPLE.value_counts()

        snv_summary = snv_counts.describe()
        indel_summary = indel_counts.describe()
        snv_q1 = snv_summary['25%']
        snv_q3 = snv_summary['75%']
        indel_q1 = indel_summary['25%']
        indel_q3 = indel_summary['75%']
        snv_iqr = snv_q3 - snv_q1
        indel_iqr = indel_q3 - indel_q1

        snv_outliers = snv_counts[(snv_counts>=(snv_q3 + 1.5*snv_iqr))].index 
        indel_outliers = indel_counts[(indel_counts>=(indel_q3 + 1.5*indel_iqr))].index 

        cohort_outliers = pd.Series(np.union1d(snv_outliers, indel_outliers))
        cohort_outliers.to_csv('outliers.txt', header=None, index=False)
        EOF

        python3 get_outliers.py ~{cohort_prefix} ~{vcf_metrics_tsv}
    >>>

    output {
        Array[String] outlier_samples = read_lines('outliers.txt')
    }
}