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

workflow filterHTtoVCF {
    input {
        String input_ht
        String hail_docker
        String genome_build
        Float gnomad_af_threshold=0.001
        Float AC_threshold=5
    }
    
    call helpers.getHailMTSize as getInputMTSize {
        input:
            mt_uri=input_ht,
            hail_docker=hail_docker
    }

    call filterHT {
        input:
        input_ht=input_ht,
        input_size=getInputMTSize.mt_size,
        hail_docker=hail_docker,
        genome_build=genome_build,
        gnomad_af_threshold=gnomad_af_threshold,
        AC_threshold=AC_threshold
    }

    output {
        File filtered_vcf_file = filterHT.filtered_vcf_file
        File filtered_vcf_idx = filterHT.filtered_vcf_idx
    }
}

task filterHT {
    input {
        String input_ht
        String hail_docker
        String genome_build
        Float input_size
        Float gnomad_af_threshold
        Float AC_threshold
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
        set -eou pipefail
        cat <<EOF > filter_ht.py
        import datetime
        import pandas as pd
        import hail as hl
        import numpy as np
        import sys
        import os

        input_ht = sys.argv[1]
        gnomad_af_threshold = float(sys.argv[2])
        AC_threshold = int(sys.argv[3])
        genome_build = sys.argv[4]
        cores = sys.argv[5]
        mem = int(np.floor(float(sys.argv[6])))

        hl.init(min_block_size=128, 
                local=f"local[*]", 
                spark_conf={
                            "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                            "spark.speculation": 'true'
                            }, 
                tmp_dir="tmp", local_tmpdir="tmp",
                            )

        ht = hl.read_table(input_ht)

        ht = ht.annotate(info=ht.info.annotate(** {field: ht.vep.worst_csq[field] for field in list(ht.vep.worst_csq)} ))

        gnomad_fields = [field for field in list(ht.vep.worst_csq) if 'gnomAD' in field]
        ht = ht.annotate(info=ht.info.annotate(gnomad_popmax_af=hl.max([hl.or_missing(hl.array(hl.set(ht.vep.transcript_consequences[gnomad_field]))[0]!='',
                                            hl.float(hl.array(hl.set(ht.vep.transcript_consequences[gnomad_field]))[0])) 
                                    for gnomad_field in gnomad_fields])))

        filt_ht = ht.filter(ht.filters.size()==0)  # PASS only
        filt_ht = filt_ht.filter((filt_ht.info.gnomad_popmax_af<=gnomad_af_threshold) | (filt_ht.info.cohort_AC<=AC_threshold))
        output_uri = os.path.basename(input_ht).split('.ht')[0] + "_gnomAD_AF_AC_PASS_only.vcf.bgz"
        hl.export_vcf(filt_ht, output_uri, tabix=True)
        EOF
        python3 filter_ht.py ~{input_ht} ~{gnomad_af_threshold} ~{AC_threshold} ~{genome_build} ~{cpu_cores} ~{memory}
    >>>

    output {
        File filtered_vcf_file = basename(input_ht, '.ht') + "_gnomAD_AF_AC_PASS_only.vcf.bgz"
        File filtered_vcf_idx = basename(input_ht, '.ht') + "_gnomAD_AF_AC_PASS_only.vcf.bgz.tbi"
    }
}