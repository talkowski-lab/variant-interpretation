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

workflow renameVCFSamples_workflow {
    input {
        File vcf_file
        File sample_map_tsv
        String sv_base_mini_docker
    }

    call renameVCFSamples {
        input:
        vcf_file=vcf_file,
        sample_map_tsv=sample_map_tsv,
        sv_base_mini_docker=sv_base_mini_docker
    }

    output {
        File renamed_vcf = renameVCFSamples.output_vcf
        File renamed_vcf_idx = renameVCFSamples.output_vcf_idx
    }
}

task renameVCFSamples {
    input {
        File vcf_file
        File sample_map_tsv
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, 'GB')
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    String output_filename = basename(vcf_file, file_ext) + '_renamed_samples.vcf.gz'
    command {
        bcftools reheader -s ~{sample_map_tsv} -o ~{output_filename} ~{vcf_file}
        tabix ~{output_filename}
    }

    output {
        File output_vcf = output_filename
        File output_vcf_idx = output_filename + '.tbi'
    }
}
