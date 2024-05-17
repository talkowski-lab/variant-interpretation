version 1.0

import "wes-denovo-helpers.wdl" as helpers
import "scatterHailMTs.wdl" as scatterHailMTs

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow exportVDStoVCF {
    input {
        File sample_file
        String input_vds
        String hail_docker
        String output_vcf_basename
        Int n_shards
    }
    
    call helpers.getHailMTSize as getInputVDSSize {
        input:
            mt_uri=input_vds,
            hail_docker=hail_docker
    }

    call scatterHailMTs.getRepartitions as getRepartitions {
        input:
        n_shards=n_shards,
        mt_uri=input_vds,
        hail_docker=hail_docker
    }

    scatter (interval in getRepartitions.partition_intervals) {
        call exportVDS {
            input:
                sample_file=sample_file,
                input_vds=input_vds,
                output_vcf_basename=output_vcf_basename,
                shard_n=interval[0],
                interval_start=interval[1],
                interval_end=interval[2],
                hail_docker=hail_docker,
                input_size=getInputVDSSize.mt_size / n_shards
        }
    }


    output {
        Array[File] vcf_shards = exportVDS.vcf_shard
        Array[File] vcf_shards_index = exportVDS.vcf_shard_idx
    }
}

task exportVDS {
    input {
        File sample_file
        String input_vds
        String output_vcf_basename
        Int shard_n
        Int interval_start
        Int interval_end
        String hail_docker
        Float input_size
        RuntimeAttr? runtime_attr_override
    }

    Float base_disk_gb = 10.0
    Float input_disk_scale = 2.0

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
    set -eou pipefail
    cat <<EOF > export_vds.py
    import pandas as pd
    import numpy as np
    import hail as hl
    import os
    import sys

    input_vds = sys.argv[1]
    output_vcf_basename = sys.argv[2]
    sample_file = sys.argv[3]
    shard_n = int(sys.argv[4])
    interval_start = int(sys.argv[5])
    interval_end = int(sys.argv[6])
    cores = sys.argv[7]  # string
    mem = int(np.floor(float(sys.argv[8])))

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{mem}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{mem}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    samples = pd.read_csv(sample_file, header=None)[0].tolist()

    vds = hl.vds.read_vds(input_vds)
    mt = vds.variant_data

    # get specific range of partitions
    mt = mt._filter_partitions(range(interval_start, interval_end))
    # subset samples
    mt = mt.filter_cols(hl.array(samples).contains(mt.s))
    # convert LGT to GT
    mt = mt.annotate_entries(GT=hl.vds.lgt_to_gt(mt.LGT, mt.LA))

    # remove all AC=0
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.AC[1] > 0, keep = True)
    mt = mt.drop('variant_qc')
    
    # move gvcf_info from entries to rows
    rows = mt.entries().select('rsid','gvcf_info').key_by('locus', 'alleles')
    mt = mt.annotate_rows(info=rows[mt.row_key].gvcf_info).drop('gvcf_info')

    hl.export_vcf(mt, f"{output_vcf_basename}_shard_{shard_n}.vcf.bgz", tabix=True)
    EOF
    python3 export_vds.py ~{input_vds} ~{output_vcf_basename} ~{sample_file} \
        ~{shard_n} ~{interval_start} ~{interval_end} ~{cpu_cores} ~{memory}
    >>>

    output {
        File vcf_shard = "~{output_vcf_basename}_shard_~{shard_n}.vcf.bgz"
        File vcf_shard_idx = "~{output_vcf_basename}_shard_~{shard_n}.vcf.bgz.tbi"
    }
}
