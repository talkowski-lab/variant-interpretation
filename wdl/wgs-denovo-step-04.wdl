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
        File ped_sex_qc
        Array[File] split_trio_vcfs
        String get_sample_pedigree_script
        String trio_denovo_docker
        Float minDQ
    }
    scatter (vcf_file in split_trio_vcfs) {
        call trio_denovo {
            input:
                ped_sex_qc=ped_sex_qc,
                vcf_file=vcf_file,
                get_sample_pedigree_script=get_sample_pedigree_script,
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

task trio_denovo {
    input {
        File ped_sex_qc
        File vcf_file
        String get_sample_pedigree_script
        String trio_denovo_docker
        Float minDQ
    }

    runtime {
        docker: trio_denovo_docker
    }

    command <<<
        sample=$(basename ~{vcf_file} '.vcf' | awk -F "_trio_" '{print $2}') 
        sample="${sample//_HP_VAF/}"
        curl ~{get_sample_pedigree_script} > get_sample_pedigree_script.py
        python3 get_sample_pedigree_script.py ~{ped_sex_qc} $sample
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
