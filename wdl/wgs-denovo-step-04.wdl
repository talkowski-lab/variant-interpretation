version 1.0 

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow step4 {
    input {
        File ped_uri
        Array[File] split_trio_vcfs
        File fill_missing_pl_script
        File get_sample_pedigree_py
        String hail_docker
        String trio_denovo_docker
        Float minDQ
        RuntimeAttr? runtime_attr_fill_missing_pl
    }
    scatter (vcf_file in split_trio_vcfs) {
        # call fillMissingPL {
        #     input:
        #         fill_missing_pl_script=fill_missing_pl_script,
        #         vcf_file=vcf_file,
        #         hail_docker=hail_docker,
        #         runtime_attr_override=runtime_attr_fill_missing_pl
        # }

        call trio_denovo {
            input:
                ped_uri=ped_uri,
                vcf_file=vcf_file,
                get_sample_pedigree_py=get_sample_pedigree_py,
                trio_denovo_docker=trio_denovo_docker,
                minDQ=minDQ
        }
    }
    call combineOutputVCFs {
        input:
            out_vcfs=trio_denovo.out_vcf,
            trio_denovo_docker=trio_denovo_docker
    }

    output {
        Array[File] trio_denovo_vcf = combineOutputVCFs.trio_denovo_vcf
    }
}

task fillMissingPL {
    input {
        File fill_missing_pl_script
        File vcf_file
        String hail_docker
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
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String out_vcf_filename = basename(vcf_file)

    command {
        python3 ~{fill_missing_pl_script} ~{vcf_file} ~{out_vcf_filename}
    }

    output {
        File out_vcf = out_vcf_filename
    }
}

task trio_denovo {
    input {
        File ped_uri
        File vcf_file
        File get_sample_pedigree_py
        String trio_denovo_docker
        Float minDQ
    }

    runtime {
        docker: trio_denovo_docker
    }

    command <<<
        sample=$(basename ~{vcf_file} '.vcf' | awk -F "_trio_" '{print $2}') 
        sample="${sample//_HP_VAF/}"
        python3 ~{get_sample_pedigree_py} ~{ped_uri} $sample
        /src/wgs_denovo/triodenovo/triodenovo-fix/src/triodenovo --ped "$sample".ped \
            --in_vcf ~{vcf_file} \
            --out_vcf ~{basename(vcf_file, '.vcf') + '.denovos.vcf'} \
            --minDQ ~{minDQ}
        bgzip ~{basename(vcf_file, '.vcf') + '.denovos.vcf'}
    >>>

    output {
        File out_vcf = basename(vcf_file, '.vcf') + '.denovos.vcf.gz'
    }
}

task combineOutputVCFs {
    input {
        Array[File] out_vcfs
        String trio_denovo_docker
    }

    runtime {
        docker: trio_denovo_docker
    }

    command {
        mkdir -p tmp_out_vcfs
        mv ~{sep=" " out_vcfs} tmp_out_vcfs/
    }

    output {
        Array[File] trio_denovo_vcf = glob('tmp_out_vcfs/*')
    }
}
