version 1.0

import "wes-denovo-helpers.wdl" as helpers
import "scatterVCF.wdl" as scatterVCF

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
        Array[String]? contigs
        File contig_lengths_file
        String split_vcf_hail_script
        String cohort_prefix
        String hail_docker
        Int records_per_shard
        Boolean localize_vcf
        RuntimeAttr? runtime_attr_split_into_shards
    }
    if (defined(chrom_shards_file)) {
        Array[File] chromosome_shards_from_file = read_lines(select_first([chrom_shards_file]))
    }
    Array[File] chromosome_shards_ = select_first([chromosome_shards, chromosome_shards_from_file])

    Map[String, Float] contig_lengths = read_map(contig_lengths_file)
    Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
    Array[String] contigs_ = select_first([contigs, chromosomes])
    Array[Pair[String, File]] split_chromosomes  = zip(contigs_, chromosome_shards_)

    scatter (chrom_pair in split_chromosomes) {
        String chromosome = chrom_pair.left
        File chrom_shard = chrom_pair.right
        Float chrom_n_records = contig_lengths[chromosome]
        Int n_shards = ceil(chrom_n_records / records_per_shard)
        if (localize_vcf) {
            call scatterVCF.scatterVCF as scatterChromosomes {
                input:
                    vcf_file=chrom_shard,
                    split_vcf_hail_script=split_vcf_hail_script,
                    n_shards=n_shards,
                    hail_docker=hail_docker,
                    runtime_attr_override=runtime_attr_split_into_shards
            }
        }
        if (!localize_vcf) {
            String mt_uri = chrom_shard
            call helpers.getHailMTSize as getHailMTSize {
                input:
                    mt_uri=mt_uri,
                    hail_docker=hail_docker
            }
            call scatterVCF.scatterVCFRemote as scatterChromosomesRemote {
                input:
                    vcf_file=mt_uri,
                    input_size=getHailMTSize.mt_size,
                    split_vcf_hail_script=split_vcf_hail_script,
                    n_shards=select_first([n_shards]),
                    records_per_shard=select_first([records_per_shard, 0]),
                    hail_docker=hail_docker,
                    runtime_attr_override=runtime_attr_split_into_shards
            }
        }
        Array[File] chrom_shards = select_first([scatterChromosomes.shards, scatterChromosomesRemote.shards])
    }

    Array[File] shards = flatten(chrom_shards)

    output {
    Array[File] vcf_shards = shards
    }
}
