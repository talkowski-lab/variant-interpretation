version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow relatedness {
    input {
        File vcf_file
        String relatedness_docker
        RuntimeAttr? runtime_attr_relatedness  
        RuntimeAttr? runtime_attr_sex  
    }

    call VCFToolsRelatedness {
        input:
            vcf_file=vcf_file,
            relatedness_docker=relatedness_docker,
            runtime_attr_override=runtime_attr_relatedness
    }

    call ImputeSexPLINK {
        input:
            vcf_file=vcf_file,
            relatedness_docker=relatedness_docker,
            runtime_attr_override=runtime_attr_sex
    }

    output {
        File out_relatedness = VCFToolsRelatedness.out_relatedness
        File out_sex = ImputeSexPLINK.out_sex
    }
}

task VCFToolsRelatedness {
    input {
        File vcf_file
        String relatedness_docker
        RuntimeAttr? runtime_attr_override  
    }

    Float relatedness_size = size(vcf_file, "GB") 
    Float base_disk_gb = 10.0
    RuntimeAttr runtime_default = object {
                                      mem_gb: 16,
                                      disk_gb: ceil(base_disk_gb + (relatedness_size) * 5.0),
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
        docker: relatedness_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    output{
        File out_relatedness = "out.relatedness2"
    }

    command <<<
        ##Run VCFtools
        vcftools --gzvcf ~{vcf_file} --relatedness2
    >>>
}

task ImputeSexPLINK {
    input {
        File vcf_file
        String relatedness_docker
        RuntimeAttr? runtime_attr_override  
    }

    Float relatedness_size = size(vcf_file, "GB") 
    Float base_disk_gb = 10.0
    RuntimeAttr runtime_default = object {
                                      mem_gb: 16,
                                      disk_gb: ceil(base_disk_gb + (relatedness_size) * 5.0),
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
        docker: relatedness_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command {
        plink --vcf ~{vcf_file} --split-x 'hg38' --make-bed --out ~{basename(vcf_file, '.vcf.gz') + '_split_x'}
        plink --bfile ~{basename(vcf_file, '.vcf.gz') + '_split_x'} --recode vcf bgz --out ~{basename(vcf_file, '.vcf.gz') + '_split_x'}
        plink --vcf ~{basename(vcf_file, '.vcf.gz') + '_split_x.vcf.gz'} --check-sex
    }

    output {
        File out_sex = 'plink.sexcheck'
    }
}