version 1.0

import "wes-denovo-helpers.wdl" as helpers
import "scatterHailMTs.wdl" as scatterHailMTs

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow exportVDStoVCF {
    input {
        File sample_file
        Array[File] vcf_shards
        String info_ht_uri
        String qc_ht_uri
        String vep_ht_uri
        String hail_docker
        String add_info_to_vcf_script
    }

    scatter (vcf_uri in vcf_shards) {
        call addInfo {
            input:
                vcf_uri=vcf_uri,
                info_ht_uri=info_ht_uri,
                qc_ht_uri=qc_ht_uri,
                vep_ht_uri=vep_ht_uri,
                hail_docker=hail_docker,
                add_info_to_vcf_script=add_info_to_vcf_script
        }
    }

    output {
        Array[File] vcf_shard_with_info = addInfo.vcf_shard
        Array[File] vcf_shards_with_info_index = addInfo.vcf_shard_idx
    }
}

task addInfo {
    input {
        File vcf_uri
        String info_ht_uri
        String qc_ht_uri
        String vep_ht_uri
        String add_info_to_vcf_script
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_uri, 'GB')
    Float base_disk_gb = 10.0
    Float input_disk_scale = 2.0

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
    set -eou pipefail
    curl ~{add_info_to_vcf_script} > add_info.py
    python3 add_info.py ~{vcf_uri} ~{info_ht_uri} ~{vep_ht_uri} ~{qc_ht_uri} \
        ~{cpu_cores} ~{memory}
    >>>

    output {
        File vcf_shard = "~{basename(vcf_uri, '.vcf.bgz')}_info.vcf.bgz"
        File vcf_shard_idx = "~{basename(vcf_uri, '.vcf.bgz')}_info.vcf.bgz.tbi"
    }
}
