version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow MergeVCFs {
    input {
        Array[File] vcf_files
        String sample_set_id
        String sv_base_mini_docker
    }

    call mergeVCFs {
        input:
        vcf_files=vcf_files,
        output_vcf_name=sample_set_id + '.merged.vcf.gz',
        sv_base_mini_docker=sv_base_mini_docker
    }

    output {
        File merged_vcf_file = mergeVCFs.merged_vcf_file
        File merged_vcf_idx = mergeVCFs.merged_vcf_idx
    }
}

task mergeVCFs {
    input {
        Array[File] vcf_files
        String output_vcf_name
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_files, 'GB')
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

    command <<<
        set -euo pipefail
        VCFS="~{write_lines(vcf_files)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        for vcf in $(cat vcfs_sorted.list);
        do
            tabix $vcf
        done
        bcftools merge -m id -Oz -o ~{output_vcf_name} --file-list vcfs_sorted.list
        tabix ~{output_vcf_name}
    >>>

    output {
        File merged_vcf_file = output_vcf_name
        File merged_vcf_idx = output_vcf_name + '.tbi'
    }
}