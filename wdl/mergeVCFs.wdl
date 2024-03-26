
version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow mergeVCFs_workflow {
    input {
        Array[File] vcf_files
        String sv_base_mini_docker
        String cohort_prefix
        Boolean sort_after_merge
    }

    call mergeVCFs {
        input:
        vcf_files=vcf_files,
        sv_base_mini_docker=sv_base_mini_docker,
        cohort_prefix=cohort_prefix,
        sort_after_merge=sort_after_merge
    }

    output {
        File merged_vcf_file = mergeVCFs.merged_vcf_file
        File merged_vcf_idx = mergeVCFs.merged_vcf_idx
    }
}

task mergeVCFs {
    input {
        Array[File] vcf_files
        String sv_base_mini_docker
        String cohort_prefix
        Boolean sort_after_merge
        RuntimeAttr? runtime_attr_override
    }

    #  generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
    #  generally assume working memory is ~3 * inputs
    #  CleanVcf5.FindRedundantMultiallelics
    Float input_size = size(vcf_files, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = if !sort_after_merge then 5.0 else 15.0
    
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
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String merged_vcf_name="~{cohort_prefix}.merged.vcf.gz"
    String sorted_vcf_name="~{cohort_prefix}.merged.sorted.vcf.gz"

    command <<<
        set -euo pipefail
        VCFS="~{write_lines(vcf_files)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        bcftools concat -n --no-version -Oz --file-list vcfs_sorted.list --output ~{merged_vcf_name}
        if [ "~{sort_after_merge}" = "true" ]; then
            mkdir -p tmp
            # bcftools sort ~{merged_vcf_name} -Oz --output ~{sorted_vcf_name} -T tmp/
            cat ~{merged_vcf_name} | zcat | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' > ~{basename(sorted_vcf_name, '.gz')}
            bgzip ~{basename(sorted_vcf_name, '.gz')}
            bcftools index -t ~{sorted_vcf_name}
        else 
            bcftools index -t ~{merged_vcf_name}
        fi

    >>>

    output {
        File merged_vcf_file = if sort_after_merge then sorted_vcf_name else merged_vcf_name
        File merged_vcf_idx = if sort_after_merge then sorted_vcf_name + ".tbi" else merged_vcf_name + ".tbi"
    }
}