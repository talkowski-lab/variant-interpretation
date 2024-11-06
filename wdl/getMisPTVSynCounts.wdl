version 1.0

import "mergeVCFs.wdl" as mergeVCFs
import "wes-denovo-helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow getMisPTVSynCounts {
    input {
        Array[File] vep_vcf_files
        File ped_uri
        String mpc_ht_uri
        String cohort_prefix
        String get_csq_counts_script
        String hail_docker
        Int ad_threshold=10
        Int gq_threshold=20
        Float ab_min=0.3
        Float ab_max=0.7
        Float af_threshold=0.01
    }

    scatter (vcf_uri in vep_vcf_files) {
        call getCSQCounts {
            input:
            vcf_uri=vcf_uri,
            ped_uri=ped_uri,
            mpc_ht_uri=mpc_ht_uri,
            ad_threshold=ad_threshold,
            gq_threshold=gq_threshold,
            ab_min=ab_min,
            ab_max=ab_max,
            af_threshold=af_threshold,
            get_csq_counts_script=get_csq_counts_script,
            hail_docker=hail_docker
        }
    }

    call combineCSQCounts {
        input:
        tsvs=getCSQCounts.mis_ptv_syn_counts,
        hail_docker=hail_docker,
        merged_filename=cohort_prefix+'_mis_syn_ptv_counts.tsv',
        input_size=size(getCSQCounts.mis_ptv_syn_counts, 'GB')
    }

    output {
        File mis_ptv_syn_counts = combineCSQCounts.merged_tsv
    }
}

task getCSQCounts {
    input {
        File vcf_uri
        File ped_uri
        String mpc_ht_uri
        Int ad_threshold
        Int gq_threshold
        Float ab_min
        Float ab_max
        Float af_threshold
        String get_csq_counts_script
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_uri, "GB")
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

    String file_ext = if sub(basename(vcf_uri), '.vcf.gz', '')!=basename(vcf_uri) then '.vcf.gz' else '.vcf.bgz'

    command {
        curl ~{get_csq_counts_script} > get_csq_counts.py
        python3 get_csq_counts.py ~{vcf_uri} ~{mpc_ht_uri} ~{ped_uri} ~{ad_threshold} ~{gq_threshold} \
        ~{ab_min} ~{ab_max} ~{cpu_cores} ~{memory} ~{af_threshold}
    }

    output {
        File mis_ptv_syn_counts = basename(vcf_uri, file_ext) + '_mis_syn_ptv_counts.tsv'
    }
}

task combineCSQCounts {
     input {
        Array[String] tsvs
        String hail_docker
        String merged_filename
        Float input_size
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

    command <<<
        cat <<EOF > merge_tsvs.py
        import pandas as pd
        import numpy as np
        import sys

        tsvs = pd.read_csv(sys.argv[1], header=None)[0].tolist()
        merged_filename = sys.argv[2]

        merged_df = pd.concat([pd.read_csv(uri, sep='\t') for uri in tsvs])
        merged_df = merged_df.groupby('SYMBOL').sum().reset_index()
        merged_df.to_csv(merged_filename, sep='\t', index=False)
        EOF

        python3 merge_tsvs.py ~{write_lines(tsvs)} ~{merged_filename} > stdout
    >>>

    output {
        File merged_tsv = merged_filename 
    }   
}