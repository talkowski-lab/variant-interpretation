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
        Array[File]? chromosome_shards
        File? chrom_shards_file
        File contig_lengths_file
        String split_vcf_hail_script
        String cohort_prefix
        String vep_hail_docker
        Int records_per_shard
        RuntimeAttr? runtime_attr_split_into_shards
    }
    if (defined(chrom_shards_file)) {
        Array[File] chromosome_shards_ = read_lines(select_first([chrom_shards_file]))
    }

    Map[String, Float] contig_lengths = read_map(contig_lengths_file)
    Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
    Array[Pair[String, File]] split_chromosomes  = zip(chromosomes, select_first([chromosome_shards, chromosome_shards_]))
    scatter (chrom_pair in split_chromosomes) {
        String chromosome = chrom_pair.left
        File chrom_shard = chrom_pair.right
        Float chrom_n_records = contig_lengths[chromosome]
        Int n_shards = ceil(chrom_n_records / records_per_shard)
        call scatterVCF as scatterChromosomes {
            input:
                vcf_file=chrom_shard,
                split_vcf_hail_script=split_vcf_hail_script,
                n_shards=n_shards,
                vep_hail_docker=vep_hail_docker,
                runtime_attr_override=runtime_attr_split_into_shards
        }
    }
    Array[File] shards = flatten(scatterChromosomes.shards)

    output {
    Array[File] vcf_shards = shards
    }
}

task scatterVCF {
    input {
        File vcf_file
        String split_vcf_hail_script
        Int n_shards
        String vep_hail_docker
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
        curl ~{split_vcf_hail_script} > split_vcf_hail.py
        python3.9 split_vcf_hail.py ~{vcf_file} ~{n_shards} ~{prefix} ~{cpu_cores} ~{memory}
        for file in $(ls ~{prefix}.vcf.bgz | grep '.bgz'); do
            shard_num=$(echo $file | cut '-' -f2);
            mv ~{prefix}.vcf.bgz/$file ~{prefix}.shard_"$shard_num".vcf.bgz
        done
    >>>
    output {
        Array[File] shards = glob("~{prefix}.shard_*.vcf.bgz")
        Array[String] shards_string = glob("~{prefix}.shard_*.vcf.bgz")
    }
}
