version 1.0
    
struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow getWGSPipelineStats {
    input {
        String BILLING_PROJECT_ID
        String WORKSPACE
        String num_variants_step03_script
        String num_variants_step04_script
        String num_variants_all_steps_script
        String terra_data_table_util_script
        String vep_hail_docker
        String hail_docker
    }

    call getCohortTSV {
        input:
            BILLING_PROJECT_ID=BILLING_PROJECT_ID,
            WORKSPACE=WORKSPACE,
            hail_docker=hail_docker,
            terra_data_table_util_script=terra_data_table_util_script
    }

    call getNumVariantsStep03 {
        input:
            cohort_tsv=getCohortTSV.cohort_tsv,
            num_variants_step03_script=num_variants_step03_script,
            vep_hail_docker=vep_hail_docker
    }

    call getNumVariantsStep04 {
        input:
            cohort_tsv=getCohortTSV.cohort_tsv,
            num_variants_step04_script=num_variants_step04_script,
            vep_hail_docker=vep_hail_docker
    }

    # call combineNumVariants {
    #     input:
    #         cohort_tsv=getCohortTSV.cohort_tsv,
    #         num_vars_step03=getNumVariantsStep03.num_vars_step03,
    #         num_vars_step04=getNumVariantsStep04.num_vars_step04,
    #         num_variants_all_steps_script=num_variants_all_steps_script,
    #         vep_hail_docker=vep_hail_docker
    # }

    output {
        String num_vars_step03_dir = sub(getNumVariantsStep03.num_vars_step03[0], 
            basename(getNumVariantsStep03.num_vars_step03[0]), "")
        String num_vars_step04_dir = sub(getNumVariantsStep04.num_vars_step04[0], 
            basename(getNumVariantsStep04.num_vars_step04[0]), "")
        # File cohort_num_variants_wgs = combineNumVariants.cohort_num_variants_wgs
    }
}

task getCohortTSV {
    input {
        String BILLING_PROJECT_ID
        String WORKSPACE
        String terra_data_table_util_script
        String hail_docker
    }

    runtime {
        docker: hail_docker
    }

    command <<<
        curl ~{terra_data_table_util_script} > terra_data_table_util.py
        echo '
        import os
        import sys
        import pandas as pd
        import numpy as np
        from terra_data_table_util import get_terra_table_to_df
        BILLING_PROJECT_ID = sys.argv[1]
        WORKSPACE = sys.argv[2]

        cohort_df = get_terra_table_to_df(project=BILLING_PROJECT_ID, workspace=WORKSPACE, table_name="cohort")
        cohort_df.to_csv("cohort_data_table.tsv", sep="\t", index=False, quotechar="'")
        ' > save_cohort_tsv.py
        python3 save_cohort_tsv.py ~{BILLING_PROJECT_ID} ~{WORKSPACE}
    >>>

    output {
        File cohort_tsv = "cohort_data_table.tsv"
    }
}

task getNumVariantsStep03 {
    input {
        File cohort_tsv
        String num_variants_step03_script
        String vep_hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(cohort_tsv, "GB") 
    Float base_disk_gb = 10.0
    Float input_disk_scale = 10.0
    RuntimeAttr runtime_default = object {
        mem_gb: 8,
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
        docker: vep_hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<< 
        set -euo pipefail
        curl ~{num_variants_step03_script} > get_num_variants_step03.sh
        sed -i "s.bcftools./opt/vep/bcftools/bcftools.g" get_num_variants_step03.sh
        bash get_num_variants_step03.sh ~{cohort_tsv}
    >>>

    output {
        Array[File] num_vars_step03 = glob("num_vars_cohorts/*_split_trio_num_vars.tsv") 
    }
}

task getNumVariantsStep04 {
    input {
        File cohort_tsv
        String num_variants_step04_script
        String vep_hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(cohort_tsv, "GB") 
    Float base_disk_gb = 10.0
    Float input_disk_scale = 10.0
    RuntimeAttr runtime_default = object {
        mem_gb: 8,
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
        docker: vep_hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<< 
        set -euo pipefail
        curl ~{num_variants_step04_script} > get_num_variants_step04.sh
        sed -i "s.bcftools./opt/vep/bcftools/bcftools.g" get_num_variants_step04.sh
        bash get_num_variants_step04.sh ~{cohort_tsv}
    >>>

    output {
        Array[File] num_vars_step04 = glob("num_vars_cohorts/*_trio_denovo_num_vars.tsv") 
    }
}

task combineNumVariants {
    input {
        File cohort_tsv
        Array[File] num_vars_step03
        Array[File] num_vars_step04
        String num_variants_all_steps_script
        String vep_hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(cohort_tsv, "GB") 
    Float base_disk_gb = 10.0
    Float input_disk_scale = 10.0
    RuntimeAttr runtime_default = object {
        mem_gb: 8,
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
        docker: vep_hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        
        num_vars_step03_dir=$(dirname ~{num_vars_step03[0]})
        num_vars_step04_dir=$(dirname ~{num_vars_step04[0]})
        
        curl ~{num_variants_all_steps_script} > get_num_variants.py
        
        mkfifo /tmp/token_fifo
        ( while true ; do curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token > /tmp/token_fifo ; done ) &
        HTS_AUTH_LOCATION=/tmp/token_fifo \
        curl -sSL https://broad.io/install-gcs-connector | python3.9
        python3.9 get_num_variants.py ~{cohort_tsv} $num_vars_step03_dir $num_vars_step04_dir
    >>>

    output {
        File cohort_num_variants_wgs = "~{basename(cohort_tsv, ".tsv")}_total_avg_num_vars_per_trio.tsv"
    }
}