version 1.0

import "mergeVCFs.wdl" as mergeVCFs
import "wes-denovo-helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? gpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow filterKnownVariants {
    input {
        Array[File] vep_vcf_files
        File known_variants_file  # each line is a variant, column named 'ID'
        String cohort_prefix
        String hail_docker
        String sv_base_mini_docker
        String genome_build='GRCh38'
    }

    scatter (vcf_file in vep_vcf_files) {
        call filterVariants {
        input:
        vcf_file=vcf_file,
        known_variants_file=known_variants_file,
        hail_docker=hail_docker,
        genome_build=genome_build
        }
    }

    call mergeVCFs.mergeVCFs as mergeVCFs {
        input:
            vcf_files=filterVariants.known_vcf,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix + '_known',
            sort_after_merge=true,
            naive=false   
    }

    output {
        File merged_known_vcf = mergeVCFs.merged_vcf_file
        File merged_known_vcf_idx = mergeVCFs.merged_vcf_idx
    }
}

task filterVariants {
    input {
        File vcf_file
        File known_variants_file
        String hail_docker
        String genome_build
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

    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    
    command <<<
    cat <<EOF > filter_vcf.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    vcf_file = sys.argv[1]
    known_file = sys.argv[2]
    genome_build = sys.argv[3]
    cores = sys.argv[4]
    mem = int(np.floor(float(sys.argv[5])))

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    mt = hl.import_vcf(vcf_file, force_bgz=vcf_file.split('.')[-1] in ['.gz', '.bgz'], 
        reference_genome=genome_build, array_elements_required=False, call_fields=[])
    known_ht = hl.import_table(known_file)
    known_ht = known_ht.annotate(locus=hl.parse_variant(known_ht.ID, genome_build).locus,
                        alleles=hl.parse_variant(known_ht.ID, genome_build).alleles).key_by('locus','alleles')

    mt = mt.semi_join_rows(known_ht)
    header = hl.get_vcf_metadata(vcf_file)
    hl.export_vcf(mt, f"{os.path.basename(vcf_file).split('.vcf')[0]}.known.vcf.bgz", metadata=header)
    EOF

    python3 filter_vcf.py ~{vcf_file} ~{known_variants_file} ~{genome_build} ~{cpu_cores} ~{memory}
    >>>

    output {
        File known_vcf = basename(vcf_file, file_ext) + '.known.vcf.bgz'
    }
}