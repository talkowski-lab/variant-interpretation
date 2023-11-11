version 1.0 

workflow step4 {
    input {
        File ped_uri
        Array[Array[File]] split_trio_vcfs
        String trio_denovo_docker
        Float minDQ
    }
    scatter (vcf_files in split_trio_vcfs) {
        scatter (vcf_file in vcf_files) {
            call trio_denovo {
                input:
                    ped_uri=ped_uri,
                    vcf_file=vcf_file,
                    trio_denovo_docker=trio_denovo_docker,
                    minDQ=minDQ
            }
        }
        call combineOutputVCFs {
            input:
                out_vcfs=trio_denovo.out_vcf,
                trio_denovo_docker=trio_denovo_docker
        }
    }

    output {
        Array[Array[File]] trio_denovo_vcf = combineOutputVCFs.trio_denovo_vcf
    }
}

task trio_denovo {
    input {
        File ped_uri
        File vcf_file
        String trio_denovo_docker
        Float minDQ
    }

    runtime {
        docker: trio_denovo_docker
    }

    command <<<
        fam=$(basename ~{vcf_file} | awk -F "_trio_" '{print $1}') 
        awk -v fam="$fam" '$1==fam' ~{ped_uri} > "$fam".ped
        /src/wgs_denovo/triodenovo/triodenovo-fix/src/triodenovo --ped "$fam".ped --in_vcf ~{vcf_file} --out_vcf ~{basename(vcf_file, '.vcf') + '.denovos.vcf'} --minDQ ~{minDQ}
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
