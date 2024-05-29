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
        File old_vcf_shard
        String input_vds
        String info_ht_uri
        String qc_ht_uri
        String vep_ht_uri
        String hail_docker
        String output_vcf_basename
        String export_vds_to_vcf_script
        Int n_shards
        Float? input_size
    }

    if (!defined(input_size)) {
        call helpers.getHailMTSize as getInputVDSSize {
            input:
                mt_uri=input_vds,
                hail_docker=hail_docker
        }
    }

    call scatterHailMTs.getRepartitions as getRepartitions {
        input:
        n_shards=n_shards,
        mt_uri=input_vds,
        hail_docker=hail_docker
    }

    Float input_size_ = select_first([input_size, getInputVDSSize.mt_size])
    scatter (interval in getRepartitions.partition_intervals) {
        call exportVDS {
            input:
                sample_file=sample_file,
                input_vds=input_vds,
                old_vcf_shard=old_vcf_shard,
                info_ht_uri=info_ht_uri,
                qc_ht_uri=qc_ht_uri,
                vep_ht_uri=vep_ht_uri,
                output_vcf_basename=output_vcf_basename,
                shard_n=interval[0],
                interval_start=interval[1],
                interval_end=interval[2],
                hail_docker=hail_docker,
                export_vds_to_vcf_script=export_vds_to_vcf_script,
                input_size=input_size_
        }
    }


    output {
        Array[File] vcf_shards = exportVDS.vcf_shard
        Array[File] vcf_shards_index = exportVDS.vcf_shard_idx
    }
}

task exportVDS {
    input {
        File sample_file
        File old_vcf_shard
        String input_vds
        String info_ht_uri
        String qc_ht_uri
        String vep_ht_uri
        String output_vcf_basename

        Int shard_n
        Int interval_start
        Int interval_end
        
        String export_vds_to_vcf_script
        String hail_docker
        Float input_size
        RuntimeAttr? runtime_attr_override
    }

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
    curl ~{export_vds_to_vcf_script} > export_vds.py
    python3 export_vds.py ~{input_vds} ~{output_vcf_basename} ~{sample_file} ~{info_ht_uri} ~{vep_ht_uri} ~{qc_ht_uri} \
        ~{shard_n} ~{interval_start} ~{interval_end} ~{cpu_cores} ~{memory} ~{old_vcf_shard}
    >>>

    output {
        File vcf_shard = "~{output_vcf_basename}_shard_~{shard_n}.vcf.bgz"
        File vcf_shard_idx = "~{output_vcf_basename}_shard_~{shard_n}.vcf.bgz.tbi"
    }
}