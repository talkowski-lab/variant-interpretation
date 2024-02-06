version 1.0

import "wgs-denovo-step-01-annot.wdl" as step1
import "wgs-denovo-step-03-annot.wdl" as step3

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow annotateAllSteps {
    input {
        Array[Array[File]]? vep_annotated_final_vcf
        Array[File]? vep_vcf_files
        Array[File] split_trio_vcfs
        Array[File] trio_denovo_vcf
        File merged_preprocessed_vcf_file
        String mpc_dir
        File mpc_chr22_file
        File loeuf_file
        File info_header
        String cohort_prefix
        String transfer_annot_script
        String annotate_vcf_script
        String hail_docker
        String sv_base_mini_docker
        Boolean bad_header=false
    }

    if (defined(vep_annotated_final_vcf)) {
        Array[File] vep_annotated_final_vcf_arr = flatten(select_first([vep_annotated_final_vcf]))
    }

    File vep_uri = select_first([vep_vcf_files, vep_annotated_final_vcf_arr])[0]
    
    call step1.saveVCFHeader as saveVCFHeader {
        input:
            vcf_uri=vep_uri,
            info_header=info_header,
            bad_header=bad_header,
            file_ext='.vcf' + sub(basename(vep_uri), '.*.vcf', ''),
            sv_base_mini_docker=sv_base_mini_docker
    }

    call step1.annotateStep01 as annotateStep01 {
        input:
            vcf_uri=merged_preprocessed_vcf_file,
            vep_uri=vep_uri,
            mpc_dir=mpc_dir,
            mpc_chr22_file=mpc_chr22_file,
            loeuf_file=loeuf_file,
            header_file=saveVCFHeader.header_file,
            file_ext='.vcf' + sub(basename(merged_preprocessed_vcf_file), '.*.vcf', ''),
            sample='false',
            annotate_vcf_script=annotate_vcf_script,
            hail_docker=hail_docker
    }

    scatter (vcf_uri in split_trio_vcfs) {
        call step3.annotateStep03 as annotateStep03 {
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

    call step3.mergeResults as mergeResults {
        input:
            split_trio_annot_tsvs=annotateStep03.split_trio_annot_tsv_,
            cohort_prefix=cohort_prefix,
            hail_docker=hail_docker
    }
 
    call labelStep04 {
        input:
            trio_denovo_vcf=trio_denovo_vcf,
            merged_split_trio_tsv=mergeResults.split_trio_annot_tsv,
            file_ext='.denovos.vcf.gz',
            cohort_prefix=cohort_prefix,
            transfer_annot_script=transfer_annot_script,
            hail_docker=hail_docker
    }

    output {
        File merged_preprocessed_vcf_file_annot = annotateStep01.merged_preprocessed_vcf_file_annot
        File split_trio_annot_tsv = mergeResults.split_trio_annot_tsv
        File trio_denovo_vcf_annot = labelStep04.trio_denovo_vcf_annot
    }
}

task labelStep04 {
    input {
        Array[File] trio_denovo_vcf
        File merged_split_trio_tsv
        String file_ext
        String cohort_prefix
        String transfer_annot_script
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([trio_denovo_vcf, merged_split_trio_tsv], "GB")
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
        trio_denovo_vcf_list="~{write_lines(trio_denovo_vcf)}"
        i=0
        for trio_uri in $(cat $trio_denovo_vcf_list); do
            if [[ "$i" -eq 0 ]]; then
                cat $trio_uri | zcat | grep '#' | tail -n 1 | awk -F '\t' -v OFS='\t' '{ $(NF+1) = "SAMPLE"; print }' > ~{cohort_prefix}_trio_denovos.tsv  # $(basename $trio_uri ~{file_ext}).header.txt
            fi
            i=$((i+1))
            sample=$(basename $trio_uri ~{file_ext} | sed 's/.*_trio_//g')
            cat $trio_uri | zcat | grep -v '#' | awk -F '\t' -v sample="$sample" -v OFS='\t' '{ $(NF+1) = sample; print }' >> ~{cohort_prefix}_trio_denovos.tsv
        done
        curl ~{transfer_annot_script} > transfer_annot.py
        python3 transfer_annot.py ~{cohort_prefix}_trio_denovos.tsv ~{merged_split_trio_tsv}
    >>>

    output {
        File trio_denovo_vcf_annot = "~{cohort_prefix}_trio_denovos_annot.tsv"
    }
}


