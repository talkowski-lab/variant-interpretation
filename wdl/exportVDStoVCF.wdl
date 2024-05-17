version 1.0

import "wes-denovo-helpers.wdl" as helpers

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
        String output_vcf_filename
        Int n_shards
    }
    
    call helpers.getHailMTSize as getInputVDSSize {
        input:
            mt_uri=input_vds,
            hail_docker=hail_docker
    }

    call exportVDS {
        input:
            sample_file=sample_file,
            input_vds=input_vds,
            output_vcf_filename=output_vcf_filename,
            n_shards=n_shards,
            hail_docker=hail_docker,
            input_size=getInputVDSSize.mt_size
    }

    output {
        Array[File] vcf_shards = exportVDS.vcf_shards
    }
}

task exportVDS {
    input {
        File sample_file
        String input_vds
        String output_vcf_filename
        String hail_docker
        Int n_shards
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

    String prefix = basename(output_vcf_filename, '.vcf.bgz')

    command <<<
    set -eou pipefail
    cat <<EOF > export_vds.py
    import pandas as pd
    import numpy as np
    import hail as hl
    import os
    import sys

    input_vds = sys.argv[1]
    output_vcf_filename = sys.argv[2]
    n_shards = int(sys.argv[3])
    sample_file = sys.argv[4]
    cores = sys.argv[5]  # string
    mem = int(np.floor(float(sys.argv[6])))

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{mem}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{mem}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    samples = pd.read_csv(sample_file, header=None)[0].tolist()

    vds = hl.vds.read_vds(input_vds, n_partitions=n_shards)
    mt = vds.variant_data
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

    hl.export_vcf(mt, output_vcf_filename, parallel='header_per_shard')
    EOF
    python3 export_vds.py ~{input_vds} ~{output_vcf_filename} ~{n_shards} ~{sample_file} ~{cpu_cores} ~{memory}

    for file in $(ls ~{output_vcf_filename} | grep '.bgz'); do
        shard_num=$(echo $file | cut -d '-' -f2);
        mv ~{output_vcf_filename}/$file ~{prefix}.shard_"$shard_num".vcf.bgz
    done
    >>>

    output {
        Array[File] vcf_shards = glob("~{prefix}.shard_*.vcf.bgz")
    }
}
