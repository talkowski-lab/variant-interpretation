version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow writeVCFtoHT {
    input {
        File vcf_file
        String output_ht_path
        String hail_docker
    }

    call writeHT {
        input:
        vcf_file=vcf_file,
        output_ht_path=output_ht_path,
        hail_docker=hail_docker
    }

    output {
        File hail_log = writeHT.hail_log
    }
}

task writeHT {
    input {
        File vcf_file
        String output_ht_path 
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_file, 'GB')
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
    set -eou pipefail
    cat <<EOF > write_ht.py
    from pyspark.sql import SparkSession
    import hail as hl
    import numpy as np
    import sys
    import ast
    import os

    vcf_file = sys.argv[1]
    output_ht_path = sys.argv[2]
    cores = sys.argv[3]  # string
    mem = int(np.floor(float(sys.argv[4])))

    hl.init()
    # hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
    #                     "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
    #                     "spark.driver.cores": cores,
    #                     "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
    #                     }, tmp_dir="tmp", local_tmpdir="tmp")

    file_ext = vcf_file.split('.')[-1]

    header_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'SpliceAI']
    ht = hl.import_table(vcf_file, comment=['#'], no_header=True, force_bgz=file_ext in ['gz', 'bgz'])\
        .rename({f"f{i}": col for i, col in enumerate(header_cols)})
    ht = ht.checkpoint(os.path.basename(vcf_file).split(file_ext)[0]+'ht')

    ht = ht.annotate(locus=hl.locus(ht.CHROM, hl.int(ht.POS), 'GRCh38'),
                    alleles=hl.array([ht.REF, ht.ALT]))

    ht = ht.annotate(SYMBOL=ht.SpliceAI.split('=')[1].split('|')[1])

    ht = ht.key_by('locus','alleles','SYMBOL')
    ht.write(output_ht_path, overwrite=True)
    EOF
    python3 write_ht.py ~{vcf_file} ~{output_ht_path} ~{cpu_cores} ~{memory}
    cp $(ls . | grep hail*.log) hail_log.txt
    >>>

    output {
        File hail_log = "hail_log.txt"
    }
}