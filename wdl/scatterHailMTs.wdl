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
        Array[String]? contigs
        File contig_lengths_file
        Int records_per_shard
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

        call getRepartitions {
            input:
                n_shards=n_shards,
                mt_uri=mt_uri,
                hail_docker=hail_docker            
        }

        scatter (interval in getRepartitions.partition_intervals) {
            call getMTPartitionInterval {
                input:
                    input_size=getHailMTSize.mt_size / n_shards,
                    shard_n=interval[0],
                    interval_start=interval[1],
                    interval_end=interval[2],
                    mt_uri=mt_uri,
                    bucket_id=bucket_id,
                    hail_docker=hail_docker
            }
        }
    }

    output {
        Array[String] mt_shards = flatten(getMTPartitionInterval.mt_shard)
    }
}

task getMTPartitionInterval {
    input {
        Float input_size
        Int shard_n
        Int interval_start
        Int interval_end
        String mt_uri
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
    cat <<EOF > get_mt_interval.py
    import pandas as pd
    import hail as hl
    import numpy as np 
    import datetime
    import sys
    import os

    mt_uri = sys.argv[1]
    shard_n = int(sys.argv[2])
    interval_start = int(sys.argv[3])
    interval_end = int(sys.argv[4])
    cores = sys.argv[5]
    mem = int(np.floor(float(sys.argv[6])))
    bucket_id = sys.argv[7]

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{mem}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{mem}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    mt = hl.read_matrix_table(mt_uri)
    mt = mt._filter_partitions(range(interval_start, interval_end))

    dir_name = f"{bucket_id}/hail/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))}/{os.path.basename(mt_uri).split('.mt')[0]}_shard_{shard_n}.mt"
    mt.write(dir_name)
    pd.Series([dir_name]).to_csv('mt_uri.txt', index=False, header=None)
    EOF

    python3 get_mt_interval.py ~{mt_uri} ~{shard_n} ~{interval_start} ~{interval_end} ~{cpu_cores} ~{memory} ~{bucket_id}
    >>>

    output {
        String mt_shard = read_lines('mt_uri.txt')[0]
    }
     
}

task getRepartitions {
    input {
        Int n_shards
        String mt_uri
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb),
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
    cat <<EOF > repartition_mt.py
    import pandas as pd
    import hail as hl
    import numpy as np 
    import datetime
    import sys
    import os

    mt_uri = sys.argv[1]
    n_shards = int(sys.argv[2])
    cores = sys.argv[3]
    mem = int(np.floor(float(sys.argv[4])))

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{mem}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{mem}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    mt = hl.read_matrix_table(mt_uri)

    tot_partitions = mt.n_partitions()
    partitions_per_chunk = int(np.ceil(tot_partitions / n_shards))
    points = np.arange(0, tot_partitions, partitions_per_chunk)
    intervals = pd.DataFrame([[points[i], points[i+1]] for i in range(len(points)-1)])

    intervals.to_csv('partition_intervals.tsv', sep='\t', header=None)
    EOF

    python3 repartition_mt.py ~{mt_uri} ~{n_shards} ~{cpu_cores} ~{memory}
    >>>

    output {
        Array[Array[Int]] partition_intervals = read_tsv('partition_intervals.tsv')
    }
}