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

workflow AncestryInferenceCohortSet {
    input {
        # Array[Array[File]] vep_vcf_files
        # Array[String] cohort_prefixes

        Array[File]? ancestry_vcf_files
        File? ancestry_vcf_file_
        File gnomad_vcf_uri
        File gnomad_rf_onnx
        File pop_labels_tsv
        String gnomad_loading_ht
        String infer_ancestry_script

        String cohort_set_id
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
        # if (!defined(ancestry_vcf_files)) {
        #     scatter (pair in zip(cohort_prefixes, vep_vcf_files)) {
        #         Array[File] vcf_files = pair.right
        #         String cohort_prefix = pair.left

        #         scatter (vcf_uri in vcf_files) {
        #             call subsetVCFgnomAD {
        #                 input:
        #                 vcf_uri=vcf_uri,
        #                 hail_docker=hail_docker,
        #                 gnomad_loading_ht=gnomad_loading_ht,
        #                 runtime_attr_override=runtime_attr_subset_vcfs
        #             }        
        #         }

        #         call mergeVCFs.mergeVCFs as mergeCohortVCFs {
        #             input:
        #             vcf_files=subsetVCFgnomAD.subset_vcf,
        #             sv_base_mini_docker=sv_base_mini_docker,
        #             cohort_prefix=cohort_prefix+'_gnomad_pca_sites',
        #             runtime_attr_override=runtime_attr_merge_vcfs
        #         }
        #     }
        # }
        
        call mergeVCFSamples {
            input:
            vcf_files=select_first([ancestry_vcf_files]),
            hail_docker=hail_docker,
            prefix=cohort_set_id,
            genome_build=genome_build,
            runtime_attr_override=runtime_attr_merge_vcfs
        }
    }

    File vcf_uri = select_first([ancestry_vcf_file_, mergeVCFSamples.merged_vcf_file])

    if (infer_ancestry) {
        call inferAncestry {
            input:
                vcf_uri=vcf_uri,
                gnomad_vcf_uri=gnomad_vcf_uri,
                gnomad_rf_onnx=gnomad_rf_onnx,
                pop_labels_tsv=pop_labels_tsv,
                cohort_prefix=cohort_set_id,
                gnomad_loading_ht=gnomad_loading_ht,
                infer_ancestry_script=infer_ancestry_script,
                hail_docker=hail_docker,
                num_pcs=num_pcs,
                min_prob=min_prob,
                use_gnomad_rf=use_gnomad_rf,
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
        File vcf_uri
        String gnomad_loading_ht
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_uri, "GB")
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

    String filename = basename(vcf_uri)
    String prefix = if (sub(filename, "\\.gz", "")!=filename) then basename(vcf_uri, ".vcf.gz") else basename(vcf_uri, ".vcf.bgz")

    command <<<
    cat <<EOF > filter_sites.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os
    
    vcf_uri = sys.argv[1]
    cores = sys.argv[2]
    mem = int(np.floor(float(sys.argv[3])))
    prefix = sys.argv[4]
    gnomad_loading_ht = sys.argv[5]

    hl.init(min_block_size=128, 
            local=f"local[*]", 
            spark_conf={
                        "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                        "spark.speculation": 'true'
                        }, 
            tmp_dir="tmp", local_tmpdir="tmp",
                        )

    mt = hl.import_vcf(vcf_uri, reference_genome='GRCh38', force_bgz=True, call_fields=[], array_elements_required=False)
    gnomad_loading_ht = hl.read_table(gnomad_loading_ht)
    mt = mt.filter_rows(hl.is_defined(gnomad_loading_ht[mt.row_key]))
    output_filename = prefix + '_gnomad_pca_sites.vcf.bgz'
    hl.export_vcf(mt, output_filename)
    EOF

    python3 filter_sites.py ~{vcf_uri} ~{cpu_cores} ~{memory} ~{prefix} ~{gnomad_loading_ht} > stdout

    >>>

    output {
        File subset_vcf = prefix + '_gnomad_pca_sites.vcf.bgz'
    }
}

task mergeVCFSamples {
    input {
        Array[File] vcf_files
        String prefix
        String genome_build
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_files, "GB")
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
    set -euo pipefail
    cat <<EOF > filter_sites.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os
    
    vcf_files = sys.argv[1].str.split(',')
    cores = sys.argv[2]
    mem = int(np.floor(float(sys.argv[3])))
    build = sys.argv[4]
    prefix = sys.argv[5]

    hl.init(min_block_size=128, 
            local=f"local[*]", 
            spark_conf={
                        "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                        "spark.speculation": 'true'
                        }, 
            tmp_dir="tmp", local_tmpdir="tmp",
                        )
    for i, vcf_uri in enumerate(vcf_files):
        if i==0:
            mt = hl.import_vcf(vcf_uri, reference_genome=build, force_bgz=True, call_fields=[], array_elements_required=False)    
        else:
            cohort_mt = hl.import_vcf(vcf_uri, reference_genome=build, force_bgz=True, call_fields=[], array_elements_required=False)    
            common_entry_fields = [x for x in list(np.intersect1d(list(mt.entry), list(cohort_mt.entry))) if x not in ['PGT','PID','MIN_DP']]
            mt = mt.select_entries(*common_entry_fields).union_cols(cohort_mt.select_entries(*common_entry_fields), row_join_type='inner')

    output_filename = prefix + '_gnomad_pca_sites.vcf.bgz'
    hl.export_vcf(mt, output_filename)
    EOF

    python3 filter_sites.py ~{sep=',' vcf_files} ~{cpu_cores} ~{memory} ~{genome_build} ~{prefix} > stdout

    >>>

    output {
        File merged_vcf_file = "~{prefix}_gnomad_pca_sites.vcf.gz"
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
        ~{pop_labels_tsv} ~{num_pcs} ~{min_prob} ~{cohort_prefix} ~{cpu_cores} ~{memory} ~{use_gnomad_rf} > stdout
    >>>

    output {
        File ancestry_tsv = cohort_prefix + '_inferred_ancestry.tsv'
        File ancestry_plot = cohort_prefix + '_inferred_ancestry.png'
    }
}