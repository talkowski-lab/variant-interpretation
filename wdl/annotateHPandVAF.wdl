version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow annotateHPandVAF {
    input {
        Array[Array[File]] split_trio_vcfs
        Array[File] merged_preprocessed_vcf_files
        File hg38_reference
        File hg38_reference_fai
        File hg38_reference_dict
        String jvarkit_docker
    }

    scatter (trio_vcfs in split_trio_vcfs) {
        scatter (trio_vcf in trio_vcfs) {
            call annotateVCF {
                input:
                    trio_vcf=trio_vcf,
                    merged_preprocessed_vcf_files=merged_preprocessed_vcf_files,
                    hg38_reference=hg38_reference,
                    hg38_reference_fai=hg38_reference_fai,
                    hg38_reference_dict=hg38_reference_dict,
                    jvarkit_docker=jvarkit_docker
            }
        }
        call combineOutputVCFs {
            input:
                out_vcfs=annotateVCF.split_trio_annot_vcfs,
                jvarkit_docker=jvarkit_docker
        }
    }

    output {
        Array[Array[File]] split_trio_annot_vcfs = combineOutputVCFs.split_trio_annot_vcfs
    }
}

task annotateVCF {
    input {
        File trio_vcf
        Array[File] merged_preprocessed_vcf_files
        File hg38_reference
        File hg38_reference_fai
        File hg38_reference_dict
        String jvarkit_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(trio_vcf, "GB")
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

    String clean_vcf = basename(trio_vcf, '.vcf')+'_clean.vcf'
    String hp_vcf = basename(trio_vcf, '.vcf')+'_HP.vcf'
    String out_vcf = basename(trio_vcf, '.vcf')+'_HP_VAF.vcf'

    command <<<
        bcftools head ~{trio_vcf} > old_header.txt
        bcftools head ~{merged_preprocessed_vcf_files[0]} > og_header.txt
        head -c -1 -q og_header.txt old_header.txt > new_header.txt
        bcftools reheader -h new_header.txt ~{trio_vcf} > ~{clean_vcf}

        java -jar /opt/jvarkit/dist/jvarkit.jar vcfpolyx -R ~{hg38_reference} -o ~{hp_vcf} ~{clean_vcf}
        bcftools +fill-tags ~{hp_vcf} -Ov -o ~{out_vcf} -- -t VAF 
    >>>

    output {
        File og_header = "og_header.txt"
        File new_header = "new_header.txt"
        File split_trio_annot_vcfs = out_vcf
    }
}

task combineOutputVCFs {
    input {
        Array[File] out_vcfs
        String jvarkit_docker
    }

    runtime {
        docker: jvarkit_docker
    }

    command {
        mkdir -p tmp_out_vcfs
        mv ~{sep=" " out_vcfs} tmp_out_vcfs/
    }

    output {
        Array[File] split_trio_annot_vcfs = glob('tmp_out_vcfs/*')
    }
}
