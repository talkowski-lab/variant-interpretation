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

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 4,
        boot_disk_gb: 4,
        preemptible_tries: 3,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File out_relatedness = "out.relatedness2"
    }

    command <<<
        ##Run VCFtools
        vcftools --gzvcf ~{vcf_file} --relatedness2
    >>>

    runtime {
        cpu_cores: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        boot_disk_gb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible_tries: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        max_retries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: relatedness_docker
    }
}

task ImputeSexPLINK {
    input {
        File vcf_file
        String relatedness_docker
        RuntimeAttr? runtime_attr_override  
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 4,
        boot_disk_gb: 4,
        preemptible_tries: 3,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu_cores: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        boot_disk_gb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible_tries: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        max_retries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: relatedness_docker
    }

    command {
        plink --vcf ~{vcf_file} --split-x 'hg38' --make-bed --out ~{basename(vcf_file, '.vcf.gz') + '_split_x'}
        plink --file ~{basename(vcf_file, '.vcf.gz') + '_split_x'} --recode vcf bgz
        plink --file ~{basename(vcf_file, '.vcf.gz') + '_split_x.vcf.gz'} --check-sex
    }

    output {
        File out_sex = 'plink.sexcheck'
    }
}