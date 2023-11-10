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
        File file
        File python_trio_sample_script
        File python_preprocess_script
        File lcr_uri
        File ped_uri
        Array[Array[File]] vep_annotated_final_vcf
        String sv_base_mini_docker
        String bucket_id
        String cohort_prefix
        String hail_docker
    }

    call makeTrioSampleFiles {
        input:
            python_trio_sample_script=python_trio_sample_script,
            ped_uri=ped_uri,
            bucket_id=bucket_id,
            cohort_prefix=cohort_prefix,
            hail_docker=hail_docker
    }

    String filename = basename(file)
    # if file is vcf.gz (just one file)
    Array[File] vcf_files = if (sub(filename, ".vcf.gz", "") != filename) then [file] else read_lines(file)
    
    File meta_uri = makeTrioSampleFiles.meta_uri
    File trio_uri = makeTrioSampleFiles.trio_uri
    
    Array[Pair[File, Array[File]]] vcf_list_paired = zip(vcf_files, vep_annotated_final_vcf)

    scatter (pair in vcf_list_paired) {
        File og_vcf_file = pair.left 
        Array[File] vcf_uri_sublist = pair.right
        scatter (vcf_uri in vcf_uri_sublist) {
            call preprocessVCF {
                input:
                    python_preprocess_script=python_preprocess_script,
                    lcr_uri=lcr_uri,
                    ped_uri=ped_uri,
                    vcf_uri=vcf_uri,
                    meta_uri=meta_uri,
                    trio_uri=trio_uri,
                    hail_docker=hail_docker
            }
        }
        call mergeVCFs {
            input:
                og_vcf_uri=og_vcf_file,
                vcf_contigs=preprocessVCF.preprocessed_vcf,
                sv_base_mini_docker=sv_base_mini_docker,
                cohort_prefix=cohort_prefix
        }
    }

    output {
        File ped_uri_no_header = bucket_id + "/resources/pedigrees/" + cohort_prefix + "_no_header.ped"
        Array[File] merged_preprocessed_vcf_files = mergeVCFs.merged_vcf_file
        Array[File] merged_preprocessed_vcf_idx = mergeVCFs.merged_vcf_idx
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
        File meta_uri = "~{bucket_id}/resources/metadata/~{cohort_prefix}_sample_list.txt"
        File trio_uri = "~{bucket_id}/resources/metadata/~{cohort_prefix}_trio_list.txt"
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
        String hail_docker
    }

    runtime {
        docker: hail_docker
    }

    command <<<
    python3 ~{python_preprocess_script} ~{lcr_uri} ~{ped_uri} ~{meta_uri} ~{trio_uri} ~{vcf_uri}
    >>>

    output {
        File preprocessed_vcf = basename(vcf_uri, '.vcf.gz') + '.preprocessed.vcf.bgz'
    }
}

task mergeVCFs {
    input {
        String og_vcf_uri
        Array[File] vcf_contigs
        String sv_base_mini_docker
        String cohort_prefix
        RuntimeAttr? runtime_attr_override
    }

    #  generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
    #  generally assume working memory is ~3 * inputs
    #  CleanVcf5.FindRedundantMultiallelics
    Float input_size = size(vcf_contigs, "GB")
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

    String merged_vcf_name="~{basename(og_vcf_uri, '.vcf.gz')}.vep.merged.vcf.gz"

    command <<<
        set -euo pipefail
        VCFS="~{write_lines(vcf_contigs)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        bcftools concat -n --no-version -Oz --file-list vcfs_sorted.list --output ~{merged_vcf_name}
        # java -jar /usr/picard/picard.jar GatherVcfs -I ~{sep=' -I ' vcf_contigs} -O ~{merged_vcf_name}
        bcftools index -t ~{merged_vcf_name}
    >>>

    output {
        File merged_vcf_file=merged_vcf_name
        File merged_vcf_idx=merged_vcf_name + ".tbi"
    }
}
