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
        Array[File] vep_vcf_files
        File bed_file
        String cohort_prefix
        String hail_docker
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_subset_vcfs
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_impute_sex
        RuntimeAttr? runtime_attr_check_relatedness
        RuntimeAttr? runtime_attr_plot_relatedness
    }

    scatter (vcf_uri in vep_vcf_files) {
        call subsetVCFgnomAD {
            input:
            vcf_uri=vcf_uri,
            hail_docker=hail_docker,
            runtime_attr_override=runtime_attr_subset_vcfs
        }
    }

    call mergeVCFs.mergeVCFs as mergeVCFs {
        input:
        vcf_files=subsetVCFgnomAD.subset_vcf,
        sv_base_mini_docker=sv_base_mini_docker,
        cohort_prefix=cohort_prefix+'_gnomad_pca_sites'
    }

    output {
        File ancestry_vcf_file = mergeVCFs.merged_vcf_file
        File ancestry_vcf_idx = mergeVCFs.merged_vcf_idx
    }
}

task subsetVCFgnomAD {
    input {
        File vcf_uri
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

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{mem}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{mem}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    mt = hl.import_vcf(vcf_uri, reference_genome='GRCh38', force_bgz=True, call_fields=[], array_elements_required=False)
    v3_loading_ht = hl.experimental.load_dataset(name='gnomad_pca_variant_loadings',
                                  version='3.1',
                                  reference_genome='GRCh38',
                                  region='us-central1',
                                  cloud='gcp')
    mt = mt.filter_rows(hl.is_defined(v3_loading_ht[mt.row_key]))
    output_filename = os.path.basename(vcf_uri).split('.vcf')[0] + '_gnomad_pca_sites.vcf.bgz'
    hl.export_vcf(mt, output_filename)
    EOF

    python3 filter_sites.py ~{vcf_uri} ~{cpu_cores} ~{memory} > stdout

    >>>

    String filename = basename(vcf_uri)
    String prefix = if (sub(filename, "\\.gz", "")!=filename) then basename(vcf_uri, ".vcf.gz") else basename(vcf_uri, ".vcf.bgz")

    output {
        File subset_vcf = prefix + '_gnomad_pca_sites.vcf.bgz'
    }
}