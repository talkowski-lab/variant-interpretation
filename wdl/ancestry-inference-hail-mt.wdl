version 1.0

import "wes-denovo-helpers.wdl" as helpers
import "mergeVCFs.wdl" as mergeVCFs

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow AncestryInference {
    input {
        Array[String]? vep_mt_uris
        File? ancestry_vcf_file_
        File gnomad_vcf_uri
        File gnomad_rf_onnx
        File pop_labels_tsv
        String gnomad_loading_ht
        String infer_ancestry_script

        String cohort_prefix
        String hail_docker
        String sv_base_mini_docker

        Int num_pcs=16  # for gnomADv3
        Float min_prob=0.75

        Boolean infer_ancestry=true
        Boolean use_gnomad_rf=false

        String genome_build='GRCh38'

        RuntimeAttr? runtime_attr_subset_vcfs
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_infer_ancestry
    }

    if ((!defined(ancestry_vcf_file_)) || (ancestry_vcf_file_ == '')) {
        scatter (mt_uri in select_first([vep_mt_uris])) {
            call helpers.getHailMTSize as getInputMTSize {
                input:
                    mt_uri=mt_uri,
                    hail_docker=hail_docker
            }
            call subsetVCFgnomAD {
                input:
                mt_uri=mt_uri,
                input_size=getInputMTSize.mt_size,
                hail_docker=hail_docker,
                gnomad_loading_ht=gnomad_loading_ht,
                genome_build=genome_build,
                runtime_attr_override=runtime_attr_subset_vcfs
            }
        }

        call mergeVCFs.mergeVCFs as mergeVCFs {
            input:
            vcf_files=subsetVCFgnomAD.subset_vcf,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix+'_gnomad_pca_sites',
            runtime_attr_override=runtime_attr_merge_vcfs
        }
    }

    File vcf_uri = select_first([ancestry_vcf_file_, mergeVCFs.merged_vcf_file])

    if (infer_ancestry) {
        call inferAncestry {
            input:
                vcf_uri=vcf_uri,
                gnomad_vcf_uri=gnomad_vcf_uri,
                gnomad_rf_onnx=gnomad_rf_onnx,
                pop_labels_tsv=pop_labels_tsv,
                cohort_prefix=cohort_prefix,
                gnomad_loading_ht=gnomad_loading_ht,
                infer_ancestry_script=infer_ancestry_script,
                hail_docker=hail_docker,
                num_pcs=num_pcs,
                min_prob=min_prob,
                use_gnomad_rf=use_gnomad_rf,
                genome_build=genome_build,
                runtime_attr_override=runtime_attr_infer_ancestry
        }
    }
    output {
        File ancestry_vcf_file = vcf_uri
        File ancestry_tsv = select_first([inferAncestry.ancestry_tsv])
        File ancestry_plot = select_first([inferAncestry.ancestry_plot])
    }
}

task subsetVCFgnomAD {
    input {
        String mt_uri
        String gnomad_loading_ht
        String genome_build
        String hail_docker
        Float input_size
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

    String filename = basename(mt_uri)
    String prefix = basename(mt_uri, '.mt')

    command <<<
    cat <<EOF > filter_sites.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os
    
    mt_uri = sys.argv[1]
    cores = sys.argv[2]
    mem = int(np.floor(float(sys.argv[3])))
    prefix = sys.argv[4]
    gnomad_loading_ht = sys.argv[5]
    build = sys.argv[6]

    hl.init(min_block_size=128, 
            local=f"local[*]", 
            spark_conf={
                        "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                        "spark.speculation": 'true'
                        }, 
            tmp_dir="tmp", local_tmpdir="tmp",
                        )

    mt = hl.read_matrix_table(mt_uri)
    gnomad_loading_ht = hl.read_table(gnomad_loading_ht)
    mt = mt.filter_rows(hl.is_defined(gnomad_loading_ht[mt.row_key]))
    output_filename = prefix + '_gnomad_pca_sites.vcf.bgz'
    hl.export_vcf(mt, output_filename)
    EOF

    python3 filter_sites.py ~{mt_uri} ~{cpu_cores} ~{memory} ~{prefix} ~{gnomad_loading_ht} ~{genome_build} > stdout

    >>>

    output {
        File subset_vcf = prefix + '_gnomad_pca_sites.vcf.bgz'
    }
}

task inferAncestry {
    input {
        File vcf_uri
        File gnomad_vcf_uri
        File gnomad_rf_onnx
        File pop_labels_tsv
        String cohort_prefix
        String gnomad_loading_ht
        String infer_ancestry_script
        String hail_docker
        Int num_pcs
        Float min_prob
        Boolean use_gnomad_rf
        String genome_build
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([vcf_uri, gnomad_vcf_uri], "GB")
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
    curl ~{infer_ancestry_script} > infer_ancestry.py
    python3 infer_ancestry.py ~{vcf_uri} ~{gnomad_vcf_uri} ~{gnomad_loading_ht} ~{gnomad_rf_onnx} \
        ~{pop_labels_tsv} ~{num_pcs} ~{min_prob} ~{cohort_prefix} ~{cpu_cores} ~{memory} ~{use_gnomad_rf} ~{genome_build} > stdout
    >>>

    output {
        File ancestry_tsv = cohort_prefix + '_inferred_ancestry.tsv'
        File ancestry_plot = cohort_prefix + '_inferred_ancestry.png'
    }
}