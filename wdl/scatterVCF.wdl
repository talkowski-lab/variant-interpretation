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
        File split_vcf_hail_script
        String cohort_prefix
        String vep_hail_docker
        String sv_base_mini_docker
        Boolean localize_vcf
        Boolean split_by_chromosome
        Boolean split_into_shards 
        Int compression_level=3
        Int? n_shards  
        Int? records_per_shard
        Int? thread_num_override
        RuntimeAttr? runtime_attr_split_by_chr
        RuntimeAttr? runtime_attr_split_into_shards
    }

    # shard the VCF (if not already sharded)
    if (split_by_chromosome) {
        if (localize_vcf) {
            call splitByChromosome {
                input:
                    vcf_file=file,
                    sv_base_mini_docker=sv_base_mini_docker,
                    thread_num_override=thread_num_override,
                    compression_level=compression_level,
                    runtime_attr_override=runtime_attr_split_by_chr
            }
        }

        if (!localize_vcf) {
            String vcf_uri = file
            call splitByChromosomeRemote {
                input:
                    vcf_file=vcf_uri,
                    sv_base_mini_docker=sv_base_mini_docker,
                    thread_num_override=thread_num_override,
                    compression_level=compression_level,
                    runtime_attr_override=select_first([runtime_attr_split_by_chr])
            }
        }
    }
    if (split_into_shards) {
        # if already split into chromosomes, shard further
        if (split_by_chromosome) {
            scatter (chrom_shard in select_first([splitByChromosome.shards])) {
                File chrom_shard_basename = basename(chrom_shard)
                Int chrom_n_variants = select_first([select_first([splitByChromosomeRemote.contig_lengths])[chrom_shard_basename], 0])
                Int no_localize_n_shards = ceil(chrom_n_variants / select_first([records_per_shard]))
                call scatterVCF as scatterChromosomes {
                    input:
                        vcf_uri=chrom_shard,
                        split_vcf_hail_script=split_vcf_hail_script,
                        n_shards=select_first([n_shards, no_localize_n_shards]),
                        vep_hail_docker=vep_hail_docker,
                        localize_vcf=localize_vcf,
                        chrom_filesizes=select_first([splitByChromosomeRemote.chrom_filesizes]),
                        runtime_attr_override=runtime_attr_split_into_shards
                }
            }
            Array[File] chromosome_shards = flatten(scatterChromosomes.shards)
        }

        if (!split_by_chromosome) {
            call scatterVCF {
                input:
                    vcf_uri=file,
                    split_vcf_hail_script=split_vcf_hail_script,
                    n_shards=select_first([n_shards]),
                    vep_hail_docker=vep_hail_docker,
                    localize_vcf=localize_vcf,
                    runtime_attr_override=runtime_attr_split_into_shards
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

task splitByChromosomeRemote { 
    input {
        String vcf_file
        String sv_base_mini_docker
        Int? thread_num_override
        Int? compression_level
        RuntimeAttr? runtime_attr_override
    }
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    String prefix = basename(vcf_file, ".vcf.gz")

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: select_first([runtime_attr_override]).disk_gb,
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    Int thread_num = select_first([thread_num_override, runtime_override.cpu_cores])
    String compression_str = if defined(compression_level) then "-Oz~{compression_level}" else "-Oz"

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        chr_string="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
        echo $chr_string | tr ',' '\n' > chr_list.txt

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        for chr in $(cat chr_list.txt); do
            bcftools view ~{vcf_file} ~{compression_str} -o ~{prefix}."$chr".vcf.gz --threads ~{thread_num} -s $chr &
        done

        # get number of records in each chr
        bcftools index -s ~{vcf_file} | cut -f1,3 > contig_lengths.txt
        awk -v OFS='\t' '$1="~{prefix}."$1".vcf.gz"' contig_lengths.txt > contig_lengths_with_filenames.txt

        # get chr file sizes in GB
        ls -s --block-size=kB | grep vcf.gz | awk '{print $2"\t"$1}' > chr_file_sizes_tmp.txt
        cat chr_file_sizes_tmp.txt | cut -f1 > filenames.txt
        awk -v OFMT='%.10f' '{$2=$2/(1024^2); print $2;}' chr_file_sizes_tmp.txt > chr_file_sizes_gb.txt
        paste -d '\t' filenames.txt chr_file_sizes_gb.txt > chr_file_sizes.txt
    >>>

    output {
        Array[File] shards = glob("~{prefix}.chr*.vcf.gz")
        Array[String] shards_string = glob("~{prefix}.chr*.vcf.gz")
        Map[String, String] chrom_filesizes = read_map('chr_file_sizes.txt')
        Map[String, String] contig_lengths = read_map('contig_lengths_with_filenames.txt')
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
        String vcf_uri
        File split_vcf_hail_script
        Int n_shards
        String vep_hail_docker
        Boolean localize_vcf
        Map[String, String]? chrom_filesizes
        RuntimeAttr? runtime_attr_override
    }
    
    File vcf_file = if localize_vcf then vcf_uri else split_vcf_hail_script
    String vcf_uri_basename = basename(vcf_uri)

    Float input_size = if localize_vcf then size(select_first([vcf_file]), "GB") else select_first([chrom_filesizes])[vcf_uri_basename]
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    String prefix = basename(vcf_uri, ".vcf.gz")

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
        docker: vep_hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
    
    command <<<
        set -euo pipefail
        if [[ "~{localize_vcf}" == "false" ]]; then
            export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
            curl -sSL https://broad.io/install-gcs-connector | python3.9
        fi;
        python3.9 ~{split_vcf_hail_script} ~{select_first([vcf_file, vcf_uri])} ~{n_shards} ~{prefix} ~{cpu_cores} ~{memory}
        for file in $(ls ~{prefix}.vcf.bgz | grep '.bgz'); do
            shard_num=$(echo $file | cut -d '-' -f2);
            mv ~{prefix}.vcf.bgz/$file ~{prefix}.shard_"$shard_num".vcf.bgz
        done
    >>>
    output {
        Array[File] shards = glob("~{prefix}.shard_*.vcf.bgz")
        Array[String] shards_string = glob("~{prefix}.shard_*.vcf.bgz")
    }
}
