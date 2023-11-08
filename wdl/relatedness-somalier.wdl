version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow runSomalier {
    input {
        File sites_uri
        File hg38_fasta
        File file
        File ped_uri
        String cohort_prefix
        String somalier_docker
        RuntimeAttr? runtime_attr_relatedness
    }
    
    String filename = basename(file)
    # if file is vcf.gz (just one file)
    Array[File] vcf_files = if (sub(filename, ".vcf.gz", "") != filename) then [file] else read_lines(file)

    scatter (vcf_uri in vcf_files) {
        call relatedness {
            input:
                sites_uri=sites_uri,
                hg38_fasta=hg38_fasta,
                vcf_uri=vcf_uri,
                ped_uri=ped_uri,
                cohort_prefix=cohort_prefix,
                somalier_docker=somalier_docker,
                runtime_attr_override=runtime_attr_relatedness
        }
    }

    output {
        Array[File] out_samples = relatedness.out_samples
        Array[File] out_pairs = relatedness.out_pairs
        Array[File] out_groups = relatedness.out_groups
        Array[File] out_html = relatedness.out_html
    }
}

task relatedness {
    input {
        File sites_uri
        File hg38_fasta
        File vcf_uri
        File ped_uri
        String cohort_prefix
        String somalier_docker
        RuntimeAttr? runtime_attr_override
    }

    Float relatedness_size = size(vcf_uri, "GB") 
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
        docker: somalier_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
    command {
        somalier extract -d extracted/ --sites ~{sites_uri} -f ~{hg38_fasta} ~{vcf_uri}
        somalier relate --infer --ped ~{ped_uri} -o ~{cohort_prefix} extracted/*.somalier
    }

    output {
        File out_samples = cohort_prefix + ".samples.tsv" # creates a .ped like file with extra QC columns
        File out_pairs = cohort_prefix + "pairs.tsv" # shows IBS for all possible sample pairs
        File out_groups = cohort_prefix + ".groups.tsv" # shows pairs of samples above a certain relatedness
        File out_html = cohort_prefix + ".html" # interactive html
    }
}