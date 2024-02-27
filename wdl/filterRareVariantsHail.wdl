version 1.0

import "mergeSplitVCF.wdl" as mergeSplitVCF
import "mergeVCFs.wdl" as mergeVCFs

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow filterRareVariantsHail {
    input {
        Array[File]? vep_vcf_files
        Array[File]? vep_annotated_final_vcf
        File lcr_uri
        File ped_uri
        File? meta_uri
        File? trio_uri
        File info_header
        File hg38_reference
        File hg38_reference_dict
        File hg38_reference_fai
        String python_trio_sample_script
        String filter_rare_variants_python_script
        String hail_docker
        String vep_hail_docker
        String sv_base_mini_docker
        String jvarkit_docker
        String cohort_prefix
        Boolean merge_split_vcf
        Boolean exclude_info_filters=false
        Boolean sort_after_merge=false
        Boolean bad_header=false
        Float AF_threshold=0.005
        Int AC_threshold=2
        Int shards_per_chunk=10
        RuntimeAttr? runtime_attr_merge_chunk
        RuntimeAttr? runtime_attr_filter_vcf
        RuntimeAttr? runtime_attr_merge_results
    }  

    if (!defined(meta_uri)) {
        call makeTrioSampleFiles {
            input:
                python_trio_sample_script=python_trio_sample_script,
                ped_uri=ped_uri,
                cohort_prefix=cohort_prefix,
                hail_docker=hail_docker
        }        
    }
    File meta_uri_ = select_first([meta_uri, makeTrioSampleFiles.meta_uri])
    File trio_uri_ = select_first([trio_uri, makeTrioSampleFiles.trio_uri])

    Array[File] vep_files = select_first([vep_vcf_files, vep_annotated_final_vcf])

    if (merge_split_vcf) {
        call mergeSplitVCF.splitFile as splitVEPFiles {
            input:
                file=write_lines(vep_files),
                shards_per_chunk=shards_per_chunk,
                cohort_prefix=cohort_prefix,
                vep_hail_docker=vep_hail_docker
        }
        scatter (chunk_file in splitVEPFiles.chunks) {        
            call mergeVCFs.mergeVCFs as mergeChunk {
                input:
                    vcf_files=read_lines(chunk_file),
                    sv_base_mini_docker=sv_base_mini_docker,
                    cohort_prefix=basename(chunk_file),
                    sort_after_merge=sort_after_merge,
                    runtime_attr_override=runtime_attr_merge_chunk
            }
            call filterRareVariants {
                input:
                    vcf_file=mergeChunk.merged_vcf_file,
                    lcr_uri=lcr_uri,
                    ped_uri=ped_uri,
                    meta_uri=meta_uri_,
                    trio_uri=trio_uri_,
                    info_header=info_header,
                    filter_rare_variants_python_script=filter_rare_variants_python_script,
                    vep_hail_docker=vep_hail_docker,
                    cohort_prefix=cohort_prefix,
                    AC_threshold=AC_threshold,
                    AF_threshold=AF_threshold,
                    bad_header=bad_header,
                    exclude_info_filters=exclude_info_filters,
                    runtime_attr_override=runtime_attr_filter_vcf
            }
        }
        call mergeResults as mergeResults_merged {
            input:
                ultra_rare_variants_tsvs=filterRareVariants.ultra_rare_variants_tsv,
                vep_hail_docker=vep_hail_docker,
                cohort_prefix=cohort_prefix,
                runtime_attr_override=runtime_attr_merge_results
        }
    }
    
    if (!merge_split_vcf) {
        scatter (vcf_file in vep_files) {
            call filterRareVariants as filterRareVariants_sharded {
                input:
                    vcf_file=vcf_file,
                    lcr_uri=lcr_uri,
                    ped_uri=ped_uri,
                    meta_uri=meta_uri_,
                    trio_uri=trio_uri_,
                    info_header=info_header,
                    filter_rare_variants_python_script=filter_rare_variants_python_script,
                    vep_hail_docker=vep_hail_docker,
                    cohort_prefix=basename(vcf_file),
                    AC_threshold=AC_threshold,
                    AF_threshold=AF_threshold,
                    bad_header=bad_header,
                    exclude_info_filters=exclude_info_filters,
                    runtime_attr_override=runtime_attr_filter_vcf
                    }
        }
        call mergeResults as mergeResults_sharded {
            input:
                ultra_rare_variants_tsvs=filterRareVariants_sharded.ultra_rare_variants_tsv,
                vep_hail_docker=vep_hail_docker,
                cohort_prefix=cohort_prefix,
                runtime_attr_override=runtime_attr_merge_results
        }
    }

    File merged_ultra_rare_variants_tsv = select_first([mergeResults_merged.ultra_rare_variants_tsv, mergeResults_sharded.ultra_rare_variants_tsv])

    call convertTSVtoVCF {
        input:
        tsv=merged_ultra_rare_variants_tsv,
        vcf_file=vep_files[0],
        vep_hail_docker=vep_hail_docker
    }

    call annotatePOLYX {
        input:
        vcf_file=convertTSVtoVCF.output_vcf,
        hg38_reference=hg38_reference,
        hg38_reference_fai=hg38_reference_fai,
        hg38_reference_dict=hg38_reference_dict,
        jvarkit_docker=jvarkit_docker
    }

    output {
        # File hail_log = filterRareVariants.hail_log
        File ultra_rare_variants_tsv = merged_ultra_rare_variants_tsv
        File polyx_vcf = annotatePOLYX.polyx_vcf
    }
}

task makeTrioSampleFiles {
    input {
        String python_trio_sample_script
        File ped_uri
        String cohort_prefix
        String hail_docker
    }

    runtime {
        docker: hail_docker
    }

    command <<<
    curl ~{python_trio_sample_script} > python_trio_sample_script.py
    python3 python_trio_sample_script.py ~{ped_uri} ~{cohort_prefix} 
    >>>
    
    output {
        File meta_uri = "~{cohort_prefix}_sample_list.txt"
        File trio_uri = "~{cohort_prefix}_trio_list.txt"
        File ped_uri_no_header = "~{cohort_prefix}_no_header.ped"
    }
}

task saveVCFHeader {
    input {
        File vcf_uri
        File info_header
        String sv_base_mini_docker
        Boolean bad_header
    }

    runtime {
        docker: sv_base_mini_docker
    }

    String header_filename = basename(vcf_uri, '.vcf.gz') + '_header.txt'

    command <<<
    bcftools head ~{vcf_uri} > ~{header_filename}
    if [[ "~{bad_header}" == "true" ]]; then
        bcftools head ~{vcf_uri} | grep -v "INFO=" > no_info_header.txt
        cat no_info_header.txt ~{info_header} | sort > ~{header_filename}
    fi
    >>>

    output {
        File header_file = header_filename
    }
}

task filterRareVariants {
    input {
        File vcf_file
        File lcr_uri
        File ped_uri
        File meta_uri
        File trio_uri
        File info_header
        String filter_rare_variants_python_script
        String vep_hail_docker
        String cohort_prefix
        Int AC_threshold
        Float AF_threshold
        Boolean exclude_info_filters
        Boolean bad_header
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
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
        docker: vep_hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String header_filename = basename(vcf_file, '.vcf.gz') + '_header.txt'

    command <<<
        /opt/vep/bcftools/bcftools head ~{vcf_file} > ~{header_filename}
        if [[ "~{bad_header}" == "true" ]]; then
            /opt/vep/bcftools/bcftools head ~{vcf_file} | grep -v "INFO=" > no_info_header.txt
            cat no_info_header.txt ~{info_header} | LC_ALL=C sort > ~{header_filename}
        fi
        curl ~{filter_rare_variants_python_script} > filter_rare_variants.py
        python3.9 filter_rare_variants.py ~{lcr_uri} ~{ped_uri} ~{meta_uri} ~{trio_uri} ~{vcf_file} \
        ~{cohort_prefix} ~{cpu_cores} ~{memory} ~{AC_threshold} ~{AF_threshold} ~{exclude_info_filters} ~{header_filename}

        cp $(ls . | grep hail*.log) hail_log.txt
    >>>

    output {
        File hail_log = "hail_log.txt"
        File ultra_rare_variants_tsv = cohort_prefix + '_ultra_rare_variants.tsv'
    }
}

task mergeResults {
    input {
        Array[File] ultra_rare_variants_tsvs
        String vep_hail_docker
        String cohort_prefix
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(ultra_rare_variants_tsvs, "GB")
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
        docker: vep_hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        head -n 1 ~{ultra_rare_variants_tsvs[0]} > "~{cohort_prefix}_ultra_rare_variants.tsv"; 
        tail -n +2 -q ~{sep=' ' ultra_rare_variants_tsvs} >> "~{cohort_prefix}_ultra_rare_variants.tsv"
    >>>

    output {
        File ultra_rare_variants_tsv = cohort_prefix + '_ultra_rare_variants.tsv'
    }
}

task convertTSVtoVCF {
    input {
        File tsv
        File vcf_file  # for header
        String vep_hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(tsv, "GB") + size(vcf_file, "GB")
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
        docker: vep_hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        cat <<EOF > tsv_to_vcf.py 
        import os
        import sys
        import pandas as pd
        import numpy as np
        import hail as hl

        tsv = sys.argv[1]
        vcf_file = sys.argv[2]
        cores = sys.argv[3]
        mem = int(np.floor(float(sys.argv[4])))

        hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

        mt = hl.import_matrix_table(tsv, row_fields={'CHROM':'str','POS':'int','REF':'str','ALT':'str'})

        mt = mt.annotate_rows(locus=hl.locus(mt.CHROM, mt.POS, 'GRCh38'),
                        alleles=hl.array([mt.REF, mt.ALT]))
        mt = mt.key_rows_by('locus','alleles')

        header = hl.get_vcf_metadata(vcf_file)
        hl.export_vcf(mt, os.path.basename(tsv).split('.tsv')[0]+'.vcf', metadata=header)
        EOF

        python3.9 tsv_to_vcf.py ~{tsv} ~{vcf_file} ~{cpu_cores} ~{memory}
    >>>

    output {
        File output_vcf = basename(tsv, '.tsv') + '.vcf'
    }
}

task annotatePOLYX {
    input {
        File vcf_file
        File hg38_reference
        File hg38_reference_fai
        File hg38_reference_dict
        String jvarkit_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 2.0
    Float input_mem_scale = 3.0
    Float input_disk_scale = 5.0
    
    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb + input_size * input_mem_scale,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: jvarkit_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String out_vcf = basename(vcf_file, '.vcf')+'_POLYX.vcf'

    command {
        java -jar /opt/jvarkit/dist/jvarkit.jar vcfpolyx -R ~{hg38_reference} -o ~{out_vcf} ~{vcf_file}
    }

    output {
        File polyx_vcf = out_vcf
    }
}