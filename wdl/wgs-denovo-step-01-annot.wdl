version 1.0

import "wgs-denovo-annotate.wdl" as annotateVCF

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow annotateStep1 {
    input {
        Array[File]? vep_annotated_final_vcf
        Array[File]? vep_vcf_files
        File merged_preprocessed_vcf_file
        String mpc_dir
        File mpc_chr22_file
        File loeuf_file
        File info_header
        String cohort_prefix
        String annotate_vcf_script
        String hail_docker
        String sv_base_mini_docker
        Boolean bad_header=false
    }

    if (defined(vep_annotated_final_vcf)) {
        Array[File] vep_annotated_final_vcf_arr = flatten(select_first([vep_annotated_final_vcf]))
    }
    File vep_uri = select_first([vep_vcf_files, vep_annotated_final_vcf_arr])[0]
    
    call saveVCFHeader {
        input:
            vcf_uri=vep_uri,
            info_header=info_header,
            bad_header=bad_header,
            file_ext='.vcf' + sub(basename(vep_uri), '.*.vcf', ''),
            sv_base_mini_docker=sv_base_mini_docker
    }

    call annotateVCF.annotateVCF as annotateStep01 {
        input:
            vcf_uri=merged_preprocessed_vcf_file,
            vep_uri=vep_uri,
            header_file=saveVCFHeader.header_file,
            mpc_dir=mpc_dir,
            mpc_chr22_file=mpc_chr22_file,
            loeuf_file=loeuf_file,
            file_ext='.vcf' + sub(basename(merged_preprocessed_vcf_file), '.*.vcf', ''),
            sample='false',
            annotate_vcf_script=annotate_vcf_script,
            hail_docker=hail_docker
    }

    output {
        File merged_preprocessed_vcf_file_annot = annotateStep01.annotated_tsv
    }
}

task saveVCFHeader {
    input {
        File vcf_uri
        File info_header
        String file_ext
        String sv_base_mini_docker
        Boolean bad_header
    }

    runtime {
        docker: sv_base_mini_docker
    }

    String header_filename = basename(vcf_uri, file_ext) + '_header.txt'

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

