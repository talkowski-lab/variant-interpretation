version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow checkKnownVariants {
    input {
        File bed_file
        Array[Array[File]] vep_annotated_final_vcf  # check VEP
        Array[File] merged_preprocessed_vcf_files  # check step1-2
        Array[File] merged_preprocessed_vcf_idx  # check step1-2
        Array[Array[File]] split_trio_vcfs  # check step3
        Boolean check_vep 
        Boolean check_step1
        Boolean check_step3
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_merge_vcfs
    }

    if (check_vep) {
        scatter (vcf_files in vep_annotated_final_vcf) {
            scatter (vcf_file in vcf_files) {
                call checkKnownVariantsVCF as checkVEP {
                    input:
                        bed_file=bed_file,
                        vcf_file=vcf_file,
                        sv_base_mini_docker=sv_base_mini_docker,
                        new_filename=basename(vcf_file, '.vcf.gz')+'.filtered.vcf.gz'
                }
            }
            call mergeVCFs as mergeVEP {
                input:
                    vcf_contigs=checkVEP.filtered_vcf,
                    sv_base_mini_docker=sv_base_mini_docker, 
                    cohort_prefix=basename(bed_file, '.bed')+'filtered.vep',
                    runtime_attr_override=runtime_attr_merge_vcfs
            }
        }
    }

    if (check_step1) {
        scatter (vcf_file in merged_preprocessed_vcf_files) {
            call checkKnownVariantsVCF as checkStep1 {
                input:
                    bed_file=bed_file,
                    vcf_file=vcf_file,
                    sv_base_mini_docker=sv_base_mini_docker,
                    new_filename=basename(vcf_file, '.vcf.gz')+'.filtered.step1.vcf.gz'
            }
        }
    }

    if (check_step3) {
        scatter (vcf_files in split_trio_vcfs) {
            scatter (vcf_file in vcf_files) {
                call checkKnownVariantsVCF as checkStep3 {
                    input:
                        bed_file=bed_file,
                        vcf_file=vcf_file,
                        sv_base_mini_docker=sv_base_mini_docker,
                        new_filename=basename(vcf_file, '.vcf')+'.filtered.step3.vcf'
                }
            }
            call mergeVCFs as mergeStep3 {
                input:
                    vcf_contigs=checkStep3.filtered_vcf,
                    sv_base_mini_docker=sv_base_mini_docker, 
                    cohort_prefix=basename(bed_file, '.bed')+'filtered.step3',
                    runtime_attr_override=runtime_attr_merge_vcfs
            }
        }
    }

    output {
        Array[File] known_variants_vep = select_first([mergeVEP.merged_vcf_file])
        Array[File] known_variants_step1 = select_first([checkStep1.filtered_vcf])
        Array[File] known_variants_step3 = select_first([mergeStep3.merged_vcf_file])
    }
}

task checkKnownVariantsVCF {
    input {
        File bed_file
        File vcf_file
        String sv_base_mini_docker
        String new_filename
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command {
        bgzip -k ~{vcf_file}
        vcf_file=$(basename vcf_file '.gz').gz
        bcftools index $vcf_file
        bcftools view -R ~{bed_file} $vcf_file -o ~{new_filename}
    }

    output {
        File filtered_vcf=new_filename
    }
}

task mergeVCFs {
    input {
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

    String merged_vcf_name="~{cohort_prefix}.merged.vcf.gz"

    command <<<
        set -euo pipefail
        VCFS="~{write_lines(vcf_contigs)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        bcftools concat -n --no-version -Oz --file-list vcfs_sorted.list --output ~{merged_vcf_name}
        bcftools index -t ~{merged_vcf_name}
    >>>

    output {
        File merged_vcf_file=merged_vcf_name
        File merged_vcf_idx=merged_vcf_name + ".tbi"
    }
}
