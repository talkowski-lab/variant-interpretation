version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow saveQC {
    input {
        File vcf
        File vcf_idx
        File bed_file
        String qc_dir
        String cohort_prefix
        String trio_denovo_docker
        RuntimeAttr? runtime_attr_index
    }

    call saveQCtsvs {
        input:
            vcf=vcf,
            vcf_idx=vcf_idx,
            bed_file=bed_file,
            cohort_prefix=cohort_prefix,
            qc_dir=qc_dir,
            trio_denovo_docker=trio_denovo_docker,
            runtime_attr_override=runtime_attr_index
    }
}

task saveQCtsvs {
    input {
        File vcf
        File vcf_idx
        File bed_file
        String cohort_prefix
        String qc_dir
        String trio_denovo_docker
        RuntimeAttr? runtime_attr_override
    }
    
    Float input_size = size(vcf, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 2.0
    Float input_mem_scale = 3.0
    Float input_disk_scale = 5.0
    
    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb + input_size * input_mem_scale,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
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
        docker: trio_denovo_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command {
        bcftools query -H -R ~{bed_file} -f '%CHROM\t%POS\t%REF\t%ALT[\t%QUAL\t%DP]\n' ~{vcf} -o ~{cohort_prefix}.qual.dp.tsv &

        bcftools query -H -R ~{bed_file} -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SOR\t%INFO/ReadPosRankSum\t%INFO/QD\t%INFO/MQ\n' ~{vcf} -o ~{cohort_prefix}.qc.tsv;

        gsutil -m cp ~{cohort_prefix}.qual.dp.tsv ~{qc_dir}
        gsutil -m cp ~{cohort_prefix}.qc.tsv ~{qc_dir}
    }
}
