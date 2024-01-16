version 1.0
    
struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow scatterVCF {

    input {
        File file
        String cohort_prefix
        String sv_base_mini_docker
        Boolean split_by_chromosome
        Boolean split_into_shards 
        Int compression_level=3
        Int? records_per_shard  
        Int? thread_num_override
        RuntimeAttr? runtime_attr_split_vcf
    }

    # shard the VCF (if not already sharded)
    if (split_by_chromosome) {
        call splitByChromosome {
            input:
                vcf_file=file,
                sv_base_mini_docker=sv_base_mini_docker,
                thread_num_override=thread_num_override,
                compression_level=compression_level,
                runtime_attr_override=runtime_attr_split_vcf
        }
    }
    if (split_into_shards) {
        # if already split into chromosomes, shard further
        if (defined(splitByChromosome.shards)) {
            scatter (chrom_shard in select_first([splitByChromosome.shards])) {
                call scatterVCF as scatterChromosomes {
                    input:
                        vcf_file=chrom_shard,
                        records_per_shard=select_first([records_per_shard]),
                        sv_base_mini_docker=sv_base_mini_docker,
                        thread_num_override=thread_num_override,
                        compression_level=compression_level,
                        runtime_attr_override=runtime_attr_split_vcf
                }
            }
            Array[File] chromosome_shards = flatten(scatterChromosomes.shards)
        }

        if (!defined(splitByChromosome.shards)) {
            call scatterVCF {
                input:
                    vcf_file=file,
                    records_per_shard=select_first([records_per_shard]),
                    sv_base_mini_docker=sv_base_mini_docker,
                    thread_num_override=thread_num_override,
                    compression_level=compression_level,
                    runtime_attr_override=runtime_attr_split_vcf
            }
        }
    }
    
    output {
    Array[File] vcf_shards = select_first([scatterVCF.shards, chromosome_shards, 
                                        splitByChromosome.shards])
    }
}   

task splitFile {
    input {
        File file
        Int shards_per_chunk
        String cohort_prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(file, "GB")
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

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        split -l ~{shards_per_chunk} ~{file} -a 4 -d "~{cohort_prefix}.shard."
    >>>

    output {
        Array[File] chunks = glob("~{cohort_prefix}.*")
    }
}

task splitByChromosome { 
    input {
        File vcf_file
        String sv_base_mini_docker
        Int? thread_num_override
        Int? compression_level
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    String prefix = basename(vcf_file, ".vcf.gz")

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    Int thread_num = select_first([thread_num_override,runtime_override.cpu_cores])
    String compression_str = if defined(compression_level) then "-Oz~{compression_level}" else "-Oz"

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        chr_string="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
        echo $chr_string | tr ',' '\n' > chr_list.txt
        awk '$0="~{prefix}."$0' chr_list.txt > filenames.txt
        paste chr_list.txt filenames.txt > chr_filenames.txt
        bcftools +scatter ~{vcf_file} ~{compression_str} -o . --threads ~{thread_num} -S chr_filenames.txt 
    >>>

    output {
        Array[File] shards = glob("~{prefix}.chr*.vcf.gz")
        Array[String] shards_string = glob("~{prefix}.chr*.vcf.gz")
    }
}

task scatterVCF {
    input {
        File vcf_file
        Int records_per_shard
        String sv_base_mini_docker
        Int? compression_level
        Int? thread_num_override
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    String prefix = basename(vcf_file, ".vcf.gz")

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    Int thread_num = select_first([thread_num_override,runtime_override.cpu_cores])
    String compression_str = if defined(compression_level) then "-Oz~{compression_level}" else "-Oz"

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
    
    command <<<
        set -euo pipefail
        #  in case the file is empty create an empty shard
        bcftools view -h ~{vcf_file} | bgzip -c > "~{prefix}.0.vcf.gz"
        bcftools +scatter ~{vcf_file} -o . ~{compression_str} -p "~{prefix}". --threads ~{thread_num} -n ~{records_per_shard}

        ls "~{prefix}".*.vcf.gz | sort -k1,1V > vcfs.list
        i=0
        while read VCF; do
          shard_no=`printf %06d $i`
          mv "$VCF" "~{prefix}.shard_${shard_no}.vcf.gz"
          i=$((i+1))
        done < vcfs.list
    >>>
    output {
        Array[File] shards = glob("~{prefix}.shard_*.vcf.gz")
        Array[String] shards_string = glob("~{prefix}.shard_*.vcf.gz")
    }
}