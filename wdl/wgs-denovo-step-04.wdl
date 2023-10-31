version 1.0 

workflow step4 {
    input {
        File ped_uri
        Array[Array[[File]] split_trio_vcfs
        String trio_denovo_docker
    }
    scatter (vcf_files in split_trio_vcfs) {
        scatter (vcf_file in vcf_files) {
            String out_vcf = basename(vcf_file, '.vcf') + '_trio_denovo.vcf'
            call trio_denovo {
                input:
                    ped_uri=ped_uri,
                    vcf_file=vcf_file,
                    trio_denovo_docker=trio_denovo_docker,
                    out_vcf=out_vcf
            }
        }
    }

    output {
        Array[Array[File]] trio_denovo_vcf = trio_denovo.trio_denovo_vcf
    }
}

task trio_denovo {
    input {
        File ped_uri
        File vcf_file
        String trio_denovo_docker
        File out_vcf
    }

    runtime {
        docker: trio_denovo_docker
    }

    command {
        /src/wgs_denovo/triodenovo/triodenovo-fix/src/triodenovo --ped ~{ped_uri} --in_vcf ~{vcf_file} --out_vcf ~{out_vcf}
    }

    output {
        File trio_denovo_vcf = out_vcf
    }
}
