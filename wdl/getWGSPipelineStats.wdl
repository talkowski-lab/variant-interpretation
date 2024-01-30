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
        File cohort_tsv
        String num_variants_step03_script
        String num_variants_step04_script
        String num_variants_all_steps_script
        String vep_hail_docker
    }

    call getNumVariantsStep03 {
        input:
            cohort_tsv=cohort_tsv,
            num_variants_step03_script=num_variants_step03_script,
            vep_hail_docker=vep_hail_docker
    }

    call getNumVariantsStep04 {
        input:
            cohort_tsv=cohort_tsv,
            num_variants_step04_script=num_variants_step04_script,
            vep_hail_docker=vep_hail_docker
    }

    call combineNumVariants {
        input:
            cohort_tsv=cohort_tsv,
            num_vars_step03=getNumVariantsStep03.num_vars_step03,
            num_vars_step04=getNumVariantsStep04.num_vars_step04,
            num_variants_all_steps_script=num_variants_all_steps_script,
            vep_hail_docker=vep_hail_docker
    }

    output {
        File cohort_num_variants_wgs = combineNumVariants.cohort_num_variants_wgs
    }
}

task getNumVariantsStep03 {
    input {
        File cohort_tsv
        String num_variants_step03_script
        String vep_hail_docker
    }

    runtime {
        docker: vep_hail_docker
    }

    command <<< 
        set -euo pipefail
        curl ~{num_variants_step03_script} > get_num_variants_step03.sh
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
    }

    runtime {
        docker: vep_hail_docker
    }

    command <<< 
        set -euo pipefail
        curl ~{num_variants_step04_script} > get_num_variants_step04.sh
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
    }

    runtime {
        docker: vep_hail_docker
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
        File cohort_num_variants_wgs = "~{basename(cohort_tsv, '.tsv')}_total_avg_num_vars_per_trio.tsv"
    }
}