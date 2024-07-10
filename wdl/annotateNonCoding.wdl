version 1.0

import "mergeVCFs.wdl" as mergeVCFs

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow annotateNonCoding {
    input {
        Array[File] vep_vcf_files
        File noncoding_bed
        String cohort_prefix
        Boolean sort_after_merge
        String sv_base_mini_docker
        String hail_docker        
    }
    
    scatter (vcf_file in vep_vcf_files) {
        call annotateFromBed {
            input:
            vcf_file=vcf_file,
            noncoding_bed=noncoding_bed,
            hail_docker=hail_docker
        }
    }

    # call mergeVCFs.mergeVCFs as mergeVCFs {
    #     input:
    #     vcf_files=annotateFromBed.noncoding_vcf,
    #     sv_base_mini_docker=sv_base_mini_docker,
    #     cohort_prefix=cohort_prefix + basename(noncoding_bed, '.bed'),
    #     sort_after_merge=sort_after_merge
    # }

    output {
        Array[File] noncoding_vcf_files = annotateFromBed.noncoding_vcf
        Array[File] noncoding_vcf_idx = annotateFromBed.noncoding_vcf_idx
        # File noncoding_vcf_file = mergeVCFs.merged_vcf_file
        # File noncoding_vcf_idx = mergeVCFs.merged_vcf_idx
    }
}

task annotateFromBed {
    input {
        File vcf_file
        String noncoding_bed 
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

    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    String output_filename = basename(vcf_file, file_ext) + '_noncoding_annot' + file_ext
   
    command <<<
    cat <<EOF > annotate_noncoding.py
    from pyspark.sql import SparkSession
    import hail as hl
    import numpy as np
    import sys
    import ast
    import os

    vcf_file = sys.argv[1]
    noncoding_bed = sys.argv[2]
    cores = sys.argv[3]  # string
    mem = int(np.floor(float(sys.argv[4])))
    output_filename = sys.argv[5]

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{mem}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{mem}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    bed = hl.import_bed(noncoding_bed, reference_genome='GRCh38', skip_invalid_intervals=True)
    mt = hl.import_vcf(vcf_file, force_bgz=True, array_elements_required=False, call_fields=[], reference_genome='GRCh38')
    mt = mt.annotate_rows(info=mt.info.annotate(PREDICTED_NONCODING=bed[mt.locus].target))

    # filter only annotated
    mt = mt.filter_rows(hl.is_defined(mt.info.PREDICTED_NONCODING))

    header = hl.get_vcf_metadata(vcf_file)
    header['info']['PREDICTED_NONCODING'] = {'Description': "Class(es) of noncoding elements disrupted by SNV/Indel.", 
                                            'Number': '.', 'Type': 'String'}
    hl.export_vcf(mt, output_filename, metadata=header, tabix=True)
    EOF
    python3 annotate_noncoding.py ~{vcf_file} ~{noncoding_bed} ~{cpu_cores} ~{memory} ~{output_filename}
    >>>

    output {
        File noncoding_vcf = output_filename
        File noncoding_vcf_idx = output_filename + '.tbi'
    }
}