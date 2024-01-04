version 1.0
    
struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow vepAnnotateHailCluster {
    input {
        # File vcf_file
        # File vep_annotate_hail_python_script
        # File hg38_fasta
        # File hg38_fasta_fai
        # File human_ancestor_fa
        # File human_ancestor_fa_fai
        # File top_level_fa
        # File gerp_conservation_scores
        # File hg38_vep_cache
        # File loeuf_data
        # File mpc_file
        # String cohort_prefix
        String vep_hail_docker
        String sv_base_mini_docker
        Boolean split_vcf=true
        Int? records_per_shard
        Int? thread_num_override
        RuntimeAttr? runtime_attr_split_vcf
        RuntimeAttr? runtime_attr_vep_annotate
    }

    call testCluster {
        input:
        vep_hail_docker=vep_hail_docker
    }

    output {
        File hail_log = testCluster.hail_log
    }
}

task testCluster {
    input {
        String vep_hail_docker
    }

    runtime {
        docker: vep_hail_docker
    }

    command <<<
        echo "import hail as hl
hl.init()" > dummy.py
        gcloud config set dataproc/region us-central1
        hailctl dataproc start veptest --vep GRCh38
        hailctl dataproc submit veptest dummy.py
        hailctl dataproc stop veptest
        cp $(ls . | grep hail*.log) hail_log.txt
    >>>

    output {
        File hail_log = "hail_log.txt"
    }
}