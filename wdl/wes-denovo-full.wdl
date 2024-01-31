version 1.0

import "wes-denovo-step-01.wdl" as step1
import "wes-denovo-step-02.wdl" as step2
import "wes-denovo-step-03.wdl" as step3

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow hailDenovoWES {
    input {
        File? vcf_file
        Array[File]? vcf_files
        File ped_uri
        File purcell5k
        File mpc_chr22_file
        File loeuf_file
        Boolean sort_after_merge=false
        String mpc_dir
        String gnomad_ht_uri
        String cohort_prefix
        String hail_annotation_script
        String hail_basic_filtering_script
        String hail_denovo_filtering_script
        String hail_docker
        String sv_base_mini_docker
    }

    if (defined(vcf_files)) {
        call mergeVCFs {
            input:
                vcf_files=select_first([vcf_files]),
                sv_base_mini_docker=sv_base_mini_docker,
                cohort_prefix=cohort_prefix,
                sort_after_merge=sort_after_merge,
                merge_or_concat='concat'
        }
    }

    File vcf_file_ = select_first([vcf_file, mergeVCFs.merged_vcf_file])

    call step1.hailAnnotate as step1 {
        input:
            vcf_file=vcf_file_,
            ped_uri=ped_uri,
            purcell5k=purcell5k,
            mpc_chr22_file=mpc_chr22_file,
            mpc_dir=mpc_dir,
            gnomad_ht_uri=gnomad_ht_uri,
            cohort_prefix=cohort_prefix,
            hail_annotation_script=hail_annotation_script,
            hail_docker=hail_docker
    }

    call step2.hailBasicFiltering as step2 {
        input:
            annot_mt=step1.annot_mt,
            ped_uri=ped_uri,
            cohort_prefix=cohort_prefix,
            hail_basic_filtering_script=hail_basic_filtering_script,
            hail_docker=hail_docker
    }

    call step3.hailDenovoFiltering as step3 {
        input:
            filtered_mt=step2.filtered_mt,
            ped_uri=ped_uri,
            cohort_prefix=cohort_prefix,
            loeuf_file=loeuf_file,
            hail_denovo_filtering_script=hail_denovo_filtering_script,
            hail_docker=hail_docker
    }

    output {
        # step 1 output
        File annot_mt = step1.annot_mt
        File sample_qc_info = step1.sample_qc_info
        File pca_score_table_5k = step1.pca_score_table_5k
        File pca_loading_table_5k = step1.pca_loading_table_5k
        # step 2 output
        File filtered_mt = step2.filtered_mt
        File post_filter_sample_qc_info = step2.post_filter_sample_qc_info
        # step 3 output
        File de_novo_results = step3.de_novo_results
        File de_novo_vep = step3.de_novo_vep
        File de_novo_ht = step3.de_novo_ht
        File tdt_mt = step3.tdt_mt
        File tdt_parent_aware_mt = step3.tdt_parent_aware_mt
    }
}

task mergeVCFs {
    input {
        Array[File] vcf_files
        String sv_base_mini_docker
        String cohort_prefix
        String merge_or_concat 
        Boolean sort_after_merge
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_files, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 10.0
    
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String merged_vcf_name="~{cohort_prefix}.merged.vcf.gz"
    String sorted_vcf_name="~{cohort_prefix}.merged.sorted.vcf.gz"

    String merge_or_concat_new = if merge_or_concat == 'concat' then 'concat -n'  else merge_or_concat
    command <<<
        set -euo pipefail
        VCFS="~{write_lines(vcf_files)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        bcftools ~{merge_or_concat_new} --no-version -Oz --file-list vcfs_sorted.list --output ~{merged_vcf_name}
        if [[ "~{sort_after_merge}" == "true" ]]; then
            cat ~{merged_vcf_name} | zcat | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' > ~{basename(sorted_vcf_name, '.gz')}
            bgzip ~{basename(sorted_vcf_name, '.gz')}
            bcftools index -t ~{sorted_vcf_name}
        else
            bcftools index -t ~{merged_vcf_name}
        fi
    >>>

    output {
        File merged_vcf_file=if sort_after_merge then sorted_vcf_name else merged_vcf_name
        File merged_vcf_idx=if sort_after_merge then sorted_vcf_name else merged_vcf_name + ".tbi"
    }
}
