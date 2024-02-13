version 1.0

import "mergeSplitVCF.wdl" as mergeSplitVCF

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
        String python_trio_sample_script
        String filter_rare_variants_python_script
        String hail_docker
        String vep_hail_docker
        String sv_base_mini_docker
        String cohort_prefix
        Boolean exclude_info_filters=false
        Boolean merge_split_vcf=false
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
            call mergeVCFs as mergeChunk {
                input:
                    vcf_files=read_lines(chunk_file),
                    sv_base_mini_docker=sv_base_mini_docker,
                    cohort_prefix=basename(chunk_file),
                    sort_after_merge=sort_after_merge,
                    merge_or_concat='concat',
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

    output {
        # File hail_log = filterRareVariants.hail_log
        File ultra_rare_variants_tsv = merged_ultra_rare_variants_tsv
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

task mergeVCFs {
    input {
        Array[File] vcf_files
        String sv_base_mini_docker
        String cohort_prefix
        Boolean sort_after_merge=false
        String merge_or_concat="concat" 
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_files, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = if !sort_after_merge then 5.0 else 15.0
    
    RuntimeAttr runtime_default = object {
        mem_gb: 4,
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String merged_vcf_name="~{cohort_prefix}.merged.vcf.gz"
    String sorted_vcf_name="~{cohort_prefix}.merged.sorted.vcf.gz"

    String merge_or_concat_new = if merge_or_concat == 'concat' then 'concat -n'  else merge_or_concat
    command <<<
        set -euo pipefail
        VCFS="~{write_lines(vcf_files)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        bcftools ~{merge_or_concat_new} --no-version -Oz --file-list vcfs_sorted.list --output ~{merged_vcf_name}
        if [ "~{sort_after_merge}" = "true" ]; then
            mkdir -p tmp
            bcftools sort ~{merged_vcf_name} -Oz --output ~{sorted_vcf_name} -T tmp/
            # cat ~{merged_vcf_name} | zcat | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' > ~{basename(sorted_vcf_name, '.gz')}
            # bgzip ~{basename(sorted_vcf_name, '.gz')}
            bcftools index -t ~{sorted_vcf_name}
        else 
            bcftools index -t ~{merged_vcf_name}
        fi
    >>>

    output {
        File merged_vcf_file = if sort_after_merge then sorted_vcf_name else merged_vcf_name
        File merged_vcf_idx = if sort_after_merge then sorted_vcf_name else merged_vcf_name + ".tbi"
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