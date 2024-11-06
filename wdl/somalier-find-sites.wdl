version 1.0

import "mergeVCFs.wdl" as mergeVCFs

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow getSomalierSites {
    input {
        Array[File] vep_vcf_files
        String sv_base_mini_docker
        String cohort_prefix
    }

    call mergeVCFs.mergeVCFs as mergeVCFs {
        input:
        vcf_files=vep_vcf_files,
        sv_base_mini_docker=sv_base_mini_docker,
        cohort_prefix=cohort_prefix
    }

    call somalierFindSites {
        input:
        vcf_uri=mergeVCFs.merged_vcf_file,
        cohort_prefix=cohort_prefix,
        sv_base_mini_docker=sv_base_mini_docker
    }

    output {
        File sites_uri = somalierFindSites.sites_uri
    }
}

task somalierFindSites {
    input {
        File vcf_uri
        String cohort_prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_uri, "GB") 
    Float base_disk_gb = 10.0
    RuntimeAttr runtime_default = object {
                                      mem_gb: 16,
                                      disk_gb: ceil(base_disk_gb + (input_size) * 5.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
    command {
        bcftools index -t ~{vcf_uri}
        /somalier_test find-sites ~{vcf_uri}
        mv sites.vcf.gz ~{cohort_prefix}.sites.vcf.gz
    }

    output {
        File sites_uri = "~{cohort_prefix}.sites.vcf.gz"
    }
}