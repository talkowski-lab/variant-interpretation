
version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow step5 {
    input {
        File de_novo_merged
        File eval_regions
        String cohort_prefix
        String final_filtering_script
        String hail_docker
        Int vqslod_cutoff_snv=-20
        Int vqslod_cutoff_indel=-2
        Int AD_alt_threshold=10
        Float af_threshold=0.005
        Boolean single_variant=true
        RuntimeAttr? runtime_attr_filter_final
        RuntimeAttr? runtime_attr_eval_regions
    }

    call finalFiltering {
        input:
        de_novo_merged=de_novo_merged,
        cohort_prefix=cohort_prefix,
        final_filtering_script=final_filtering_script,
        hail_docker=hail_docker,
        vqslod_cutoff_snv=vqslod_cutoff_snv,
        vqslod_cutoff_indel=vqslod_cutoff_indel,
        af_threshold=af_threshold,
        AD_alt_threshold=AD_alt_threshold,
        single_variant=single_variant,
        runtime_attr_override=runtime_attr_filter_final
    }

    call annotateEvalRegions {
        input:
        de_novo_final=finalFiltering.de_novo_filtered,
        eval_regions=eval_regions,
        cohort_prefix=cohort_prefix,
        hail_docker=hail_docker,
        runtime_attr_override=runtime_attr_eval_regions
    }

    output {
        File de_novo_final = annotateEvalRegions.de_novo_final_eval_reg
    }
}

task finalFiltering {
    input {
        File de_novo_merged
        String cohort_prefix
        String final_filtering_script
        String hail_docker
        Int vqslod_cutoff_snv
        Int vqslod_cutoff_indel
        Int AD_alt_threshold
        Float af_threshold
        Boolean single_variant
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(de_novo_merged, 'GB')
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
    curl ~{final_filtering_script} > final_filtering.py
    python3 final_filtering.py ~{de_novo_merged} ~{cohort_prefix} ~{vqslod_cutoff_snv} ~{vqslod_cutoff_indel} \
        ~{af_threshold} ~{AD_alt_threshold} ~{cpu_cores} ~{memory} ~{single_variant} > stdout
    >>>

    output {
        File de_novo_filtered = cohort_prefix + '_de_novo_filtered_final.tsv'
    }
}

task annotateEvalRegions {
    input {
        File de_novo_final
        File eval_regions
        String cohort_prefix
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(de_novo_final, 'GB')
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
    cat <<EOF > annotate_eval_regions.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    de_novo_final = sys.argv[1]
    eval_regions = sys.argv[2]
    cores = sys.argv[3]
    mem = int(np.floor(float(sys.argv[4])))
    cohort_prefix = sys.argv[5]

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{mem}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{mem}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")
    
    ht = hl.import_table(de_novo_final)
    ht = ht.annotate(locus=hl.parse_locus(ht.locus, 'GRCh38'))
    ht = ht.key_by('locus','alleles')

    evaluation_regions = hl.import_locus_intervals(eval_regions, reference_genome = 'GRCh38')
    ht = ht.annotate(eval_reg = hl.is_defined(evaluation_regions[ht.locus]))
    
    df = ht.to_pandas()
    df.to_csv(cohort_prefix+'_eval_reg.tsv', sep='\t', index=False)
    EOF

    python3 annotate_eval_regions.py ~{de_novo_final} ~{eval_regions} ~{cpu_cores} ~{memory} ~{cohort_prefix} > stdout
    >>>

    output {
        File de_novo_final_eval_reg = cohort_prefix + '_eval_reg.tsv'
    }
}