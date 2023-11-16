version 1.0

workflow indexVCF {
    input {
        File vcf
        String vcf_dir
        String sv_base_mini_docker
    }

    call indexVCF_ {
        input:
            vcf=vcf,
            vcf_dir=vcf_dir,
            sv_base_mini_docker=sv_base_mini_docker
    }
}

task indexVCF_ {
    input {
        File vcf
        String vcf_dir
        String sv_base_mini_docker
    }
    
    runtime {
        docker: sv_base_mini_docker
    }

    command {
        tabix ~{vcf}
        gsutil -m cp ~{vcf} ~{vcf_dir}/
        gsutil -m cp ~{vcf}.tbi ~{vcf_dir}/
    }
}
