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

workflow annotateStep3 {
    input {
        Array[File]? vep_annotated_final_vcf
        Array[File]? vep_vcf_files
        Array[File] split_trio_vcfs
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

    scatter (vcf_uri in split_trio_vcfs) {
        call annotateVCF.annotateVCF as annotateStep03 {
            input:
                vcf_uri=vcf_uri,
                vep_uri=vep_uri,
                mpc_dir=mpc_dir,
                mpc_chr22_file=mpc_chr22_file,
                loeuf_file=loeuf_file,
                header_file=saveVCFHeader.header_file,
                file_ext='.vcf',
                sample=sub(basename(vcf_uri, '.vcf'), '.*_trio_', ''),
                annotate_vcf_script=annotate_vcf_script,
                hail_docker=hail_docker
        }
    }

    call mergeResults {
        input:
            annot_tsvs=annotateStep03.annotated_tsv,
            cohort_prefix=cohort_prefix,
            hail_docker=hail_docker
    }

    output {
        File split_trio_annot_tsv = mergeResults.merged_tsv
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

task mergeResults {
    input {
        Array[File] annot_tsvs
        String hail_docker
        String cohort_prefix
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(annot_tsvs, "GB")
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
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        head -n 1 ~{annot_tsvs[0]} > "~{cohort_prefix}_split_trio_vcfs_annot.tsv"; 
        tail -n +2 -q ~{sep=' ' annot_tsvs} >> "~{cohort_prefix}_split_trio_vcfs_annot.tsv"
    >>>

    output {
        File merged_tsv = cohort_prefix + '_split_trio_vcfs_annot.tsv'
    }
}