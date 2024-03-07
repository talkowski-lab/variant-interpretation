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

workflow scatterMT {
    input {
        Array[String] mt_uris
        Array[String] contigs
        File contig_lengths_file
        Int records_per_shard
        String split_mt_hail_script
        String bucket_id
        String hail_docker
    }

    Map[String, Float] contig_lengths = read_map(contig_lengths_file)
    Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
    Array[String] contigs_ = select_first([contigs, chromosomes])
    Array[Pair[String, File]] split_contigs  = zip(contigs_, mt_uris)

    scatter (contig_pair in split_contigs) {
        String contig = contig_pair.left
        String mt_uri = contig_pair.right
        Float chrom_n_records = contig_lengths[contig]
        Int n_shards = ceil(chrom_n_records / records_per_shard)

        call helpers.getHailMTSize as getHailMTSize {
            input:
                mt_uri=mt_uri,
                hail_docker=hail_docker
        }

        call scatterHailMT {
            input:
                input_size=getHailMTSize.mt_size,
                n_shards=n_shards,
                records_per_shard=records_per_shard,
                mt_uri=mt_uri,
                split_mt_hail_script=split_mt_hail_script,
                bucket_id=bucket_id,
                hail_docker=hail_docker
        }
    }

    output {
        Array[String] mt_shards = flatten(scatterHailMT.mt_shards)
    }
}

task scatterHailMT {
    input {
        Float input_size
        Int n_shards
        Int records_per_shard
        String mt_uri
        String split_mt_hail_script
        String bucket_id
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

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

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        curl ~{split_mt_hail_script} > split_mt.py
        python3 split_mt.py ~{mt_uri} ~{n_shards} ~{records_per_shard} ~{cpu_cores} ~{memory} ~{bucket_id}
    >>>

    output {
        Array[String] mt_shards = read_lines('mt_uris.txt')
    }
}
