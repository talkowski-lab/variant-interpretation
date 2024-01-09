version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow step1 {
    input {
        # file can be a list of vcf files or just one vcf file
        File python_trio_sample_script
        File python_preprocess_script
        File lcr_uri
        File ped_uri
        File info_header
        Array[Array[File]] vep_annotated_final_vcf
        Array[Array[File]] vep_annotated_final_vcf_idx
        String hail_docker
        String vep_hail_docker
        String sv_base_mini_docker
        String bucket_id
        String cohort_prefix
        Int shards_per_chunk
        Boolean bad_header=false
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_preprocess
    }

    call makeTrioSampleFiles {
        input:
            python_trio_sample_script=python_trio_sample_script,
            ped_uri=ped_uri,
            bucket_id=bucket_id,
            cohort_prefix=cohort_prefix,
            hail_docker=hail_docker
    }

    Array[File] vep_annotated_final_vcf_array = flatten(vep_annotated_final_vcf)
    Array[File] vep_annotated_final_vcf_idx_array = flatten(vep_annotated_final_vcf_idx)

    if (shards_per_chunk!=0) {
        call splitFile as splitVEPFiles {
            input:
                file=write_lines(vep_annotated_final_vcf_array),
                shards_per_chunk=shards_per_chunk,
                cohort_prefix=cohort_prefix,
                vep_hail_docker=vep_hail_docker
        }

        call splitFile as splitVEPIndexFiles {
            input:
                file=write_lines(vep_annotated_final_vcf_idx_array),
                shards_per_chunk=shards_per_chunk,
                cohort_prefix=cohort_prefix,
                vep_hail_docker=vep_hail_docker
        }

        Array[Pair[File, File]] chunk_vcf_idx_pairs = zip(splitVEPFiles.chunks, splitVEPIndexFiles.chunks)
        
        scatter (chunk_vcf_idx_pair in chunk_vcf_idx_pairs) {        
            File chunk_file = chunk_vcf_idx_pair.left    
            File chunk_idx_file = chunk_vcf_idx_pair.right 
            call mergeVCFs as mergeChunk {
                input:
                    vcf_files=read_lines(chunk_file),
                    vcf_files_idx=read_lines(chunk_idx_file),
                    sv_base_mini_docker=sv_base_mini_docker,
                    cohort_prefix=basename(chunk_file),
                    runtime_attr_override=runtime_attr_merge_vcfs
            }
            call saveVCFHeader as saveVCFHeaderChunk {
                input:
                    vcf_uri=mergeChunk.merged_vcf_file,
                    info_header=info_header,
                    bad_header=bad_header,
                    sv_base_mini_docker=sv_base_mini_docker
            }
            call preprocessVCF as preprocessVCFChunk {
                input:
                    python_preprocess_script=python_preprocess_script,
                    lcr_uri=lcr_uri,
                    ped_uri=ped_uri,
                    vcf_uri=mergeChunk.merged_vcf_file,
                    meta_uri=makeTrioSampleFiles.meta_uri,
                    trio_uri=makeTrioSampleFiles.trio_uri,
                    vep_hail_docker=vep_hail_docker,
                    header_file=saveVCFHeaderChunk.header_file,
                    runtime_attr_override=runtime_attr_preprocess
            }
        }
        call mergeVCFs as mergeChunks {
            input:
                vcf_files=preprocessVCFChunk.preprocessed_vcf,
                vcf_files_idx=preprocessVCFChunk.preprocessed_vcf_idx,
                sv_base_mini_docker=sv_base_mini_docker,
                cohort_prefix=cohort_prefix,
                runtime_attr_override=runtime_attr_merge_vcfs  
        }
    }

    if (shards_per_chunk==0) {
        scatter (vcf_uri in vep_annotated_final_vcf_array) {
            call saveVCFHeader {
                input:
                    vcf_uri=vcf_uri,
                    info_header=info_header,
                    bad_header=bad_header,
                    sv_base_mini_docker=sv_base_mini_docker
            }
            call preprocessVCF {
                input:
                    python_preprocess_script=python_preprocess_script,
                    lcr_uri=lcr_uri,
                    ped_uri=ped_uri,
                    vcf_uri=vcf_uri,
                    meta_uri=makeTrioSampleFiles.meta_uri,
                    trio_uri=makeTrioSampleFiles.trio_uri,
                    vep_hail_docker=vep_hail_docker,
                    header_file=saveVCFHeader.header_file,
                    runtime_attr_override=runtime_attr_preprocess
            }
        }
        call mergeVCFs {
            input:
                vcf_files=preprocessVCF.preprocessed_vcf,
                vcf_files_idx=preprocessVCF.preprocessed_vcf_idx,
                sv_base_mini_docker=sv_base_mini_docker,
                cohort_prefix=cohort_prefix,
                runtime_attr_override=runtime_attr_merge_vcfs
        }
    }

    output {
        File meta_uri = makeTrioSampleFiles.meta_uri
        File trio_uri = makeTrioSampleFiles.trio_uri
        File ped_uri_no_header = makeTrioSampleFiles.ped_uri_no_header
        File merged_preprocessed_vcf_file = select_first([mergeChunks.merged_vcf_file, mergeVCFs.merged_vcf_file])
        File merged_preprocessed_vcf_idx = select_first([mergeChunks.merged_vcf_idx, mergeVCFs.merged_vcf_idx])
    }
}

task makeTrioSampleFiles {
    input {
        File python_trio_sample_script
        File ped_uri
        String bucket_id
        String cohort_prefix
        String hail_docker
    }

    runtime {
        docker: hail_docker
    }

    command <<<
    python3 ~{python_trio_sample_script} ~{ped_uri} ~{cohort_prefix} ~{bucket_id}
    >>>
    
    output {
        String meta_uri = "~{bucket_id}/resources/metadata/~{cohort_prefix}_sample_list.txt"
        String trio_uri = "~{bucket_id}/resources/metadata/~{cohort_prefix}_trio_list.txt"
        String ped_uri_no_header = bucket_id + "/resources/pedigrees/" + cohort_prefix + "_no_header.ped"
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

task preprocessVCF {
    input {
        File python_preprocess_script
        File lcr_uri
        File ped_uri
        File vcf_uri
        File meta_uri
        File trio_uri
        File header_file
        String vep_hail_docker
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
        docker: vep_hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    File preprocessed_vcf_out = basename(vcf_uri, '.vcf.gz') + '.preprocessed.vcf.bgz'
    command <<<
        python3.9 ~{python_preprocess_script} ~{lcr_uri} ~{ped_uri} ~{meta_uri} ~{trio_uri} ~{vcf_uri} ~{header_file} ~{cpu_cores} ~{memory}
        /opt/vep/bcftools/bcftools index -t ~{preprocessed_vcf_out}
    >>>

    output {
        File preprocessed_vcf = preprocessed_vcf_out
        File preprocessed_vcf_idx = preprocessed_vcf_out + '.tbi'
    }
}

task mergeVCFs {
    input {
        Array[File] vcf_files
        Array[File] vcf_files_idx
        String sv_base_mini_docker
        String cohort_prefix
        RuntimeAttr? runtime_attr_override
    }

    #  generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
    #  generally assume working memory is ~3 * inputs
    #  CleanVcf5.FindRedundantMultiallelics
    Float input_size = size(vcf_files, "GB")
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String merged_vcf_name="~{cohort_prefix}.vep.merged.vcf.gz"
    String sorted_vcf_name="~{cohort_prefix}.merged.sorted.vcf.gz"

    command <<<
        set -euo pipefail
        VCFS="~{write_lines(vcf_files)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        bcftools concat -n --no-version -Oz --file-list vcfs_sorted.list --output ~{merged_vcf_name}
        bcftools sort ~{merged_vcf_name} --output ~{sorted_vcf_name}
        bcftools index -t ~{merged_vcf_name}
    >>>

    output {
        File merged_vcf_file=sorted_vcf_name
        File merged_vcf_idx=sorted_vcf_name + ".tbi"
    }
}

task splitFile {
    input {
        File file
        Int shards_per_chunk
        String cohort_prefix
        String vep_hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(file, "GB")
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

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: vep_hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        split -l ~{shards_per_chunk} ~{file} -a 4 -d "~{cohort_prefix}.shard."
    >>>

    output {
        Array[File] chunks = glob("~{cohort_prefix}.*")
    }
}
