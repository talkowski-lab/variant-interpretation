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

workflow prioritizeCSQ {
    input {
        String input_ht
        File vep_vcf_file
        String bucket_id
        String prioritize_csq_ht_script
        String hail_docker
        String sample_column
        String genome_build
    }
    
    call helpers.getHailMTSize as getInputMTSize {
        input:
            mt_uri=input_ht,
            hail_docker=hail_docker
    }

    call annotateMostSevereCSQ {
        input:
        input_ht=input_ht,
        bucket_id=bucket_id,
        vep_vcf_file=vep_vcf_file,
        prioritize_csq_ht_script=prioritize_csq_ht_script,
        hail_docker=hail_docker,
        sample_column=sample_column,
        genome_build=genome_build
    }

    output {
        File output_prior_csq = annotateMostSevereCSQ.output_prior_csq
    }
}

task annotateMostSevereCSQ {
    input {
        String input_ht
        File vep_vcf_file
        String bucket_id
        String prioritize_csq_ht_script
        String hail_docker
        String sample_column
        String genome_build
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
    
    String file_ext = if sub(basename(input_ht), '\\.gz', '')==basename(input_ht) then '.tsv' else '.tsv.gz'
    command <<<
        set -eou pipefail
        curl ~{prioritize_csq_ht_script} > prioritize_csq.py
        python3 prioritize_csq.py ~{input_ht} ~{cpu_cores} ~{memory} \
        ~{sample_column} ~{vep_vcf_file} ~{genome_build} ~{bucket_id}
    >>>

    output {
        File output_prior_csq = read_lines('ht_uri.txt')[0]
    }
}