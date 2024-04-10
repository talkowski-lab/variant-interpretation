version 1.0

import "mergeVCFs.wdl" as mergeVCFs
import "wes-denovo-helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow Relatedness {
    input {
        Array[File] vep_vcf_files
        File? merged_vep_file
        File ped_uri
        File bed_file
        Int chunk_size=0
        String cohort_prefix
        String relatedness_qc_script
        String plot_relatedness_script
        String sex_qc_script
        String sv_base_mini_docker
        String hail_docker
        String bucket_id
    }

    if (!defined(merged_vep_file)) {
        
        scatter (vcf_uri in vep_vcf_files) {
            String filename = basename(vcf_uri)
            String prefix = if (sub(filename, ".gz", "")!=filename) then basename(filename, ".vcf.gz") else basename(filename, ".vcf.bgz")
            call helpers.subsetVCFs as subsetVCFs {
                input:
                    bed_file=bed_file,
                    vcf_uri=vcf_uri,
                    vcf_idx=vcf_uri+'.tbi',
                    output_name=prefix + '.somalier.subset.vcf.gz',
                    sv_base_mini_docker=sv_base_mini_docker
            }
        }

        call mergeVCFs.mergeVCFs as mergeVCFs {
        input:
            vcf_files=subsetVCFs.subset_vcf,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix
        }
    }
    File merged_vcf_file = select_first([merged_vep_file, mergeVCFs.merged_vcf_file])

    call imputeSex {
        input:
        vcf_uri=merged_vcf_file,
        bed_file=bed_file,
        ped_uri=ped_uri,
        cohort_prefix=cohort_prefix,
        sex_qc_script=sex_qc_script,
        hail_docker=hail_docker
    }

    call checkRelatedness {
        input:
        vcf_uri=merged_vcf_file,
        bed_file=bed_file,
        ped_uri=ped_uri,
        cohort_prefix=cohort_prefix,
        relatedness_qc_script=relatedness_qc_script,
        hail_docker=hail_docker,
        bucket_id=bucket_id
    }

    call plotRelatedness {
        input:
        kinship_tsv=checkRelatedness.kinship_tsv,
        ped_uri=ped_uri,
        cohort_prefix=cohort_prefix,
        plot_relatedness_script=plot_relatedness_script,
        hail_docker=hail_docker,
        chunk_size=chunk_size
    }

    output {
        File somalier_vcf_file = merged_vcf_file
        File sex_qc_plots = imputeSex.sex_qc_plots
        File ped_sex_qc = imputeSex.ped_sex_qc
        File relatedness_qc = checkRelatedness.relatedness_qc
        File kinship_tsv = checkRelatedness.kinship_tsv
        File relatedness_plot = plotRelatedness.relatedness_plot
    }
}

task imputeSex {
    input {
        File vcf_uri
        File ped_uri
        File bed_file
        String cohort_prefix
        String sex_qc_script
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_uri, "GB")
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
        curl ~{sex_qc_script} > impute_sex.py
        python3 impute_sex.py ~{vcf_uri} ~{bed_file} ~{cohort_prefix} ~{ped_uri} ~{cpu_cores} ~{memory}
    >>>

    output {
        File sex_qc_plots = cohort_prefix + "_sex_qc.png"
        File ped_sex_qc = cohort_prefix + "_sex_qc.ped"
    }
}

task checkRelatedness {
    input {
        File vcf_uri
        File ped_uri
        File bed_file
        String cohort_prefix
        String relatedness_qc_script
        String hail_docker
        String bucket_id
        String score_table=false
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_uri, "GB")
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
        set -eou pipefail
        curl ~{relatedness_qc_script} > check_relatedness.py
        python3 check_relatedness.py ~{vcf_uri} ~{bed_file} ~{cohort_prefix} ~{ped_uri} ~{cpu_cores} ~{memory} ~{bucket_id} ~{score_table} > stdout
    >>>

    output {
        File relatedness_qc = cohort_prefix + "_relatedness_qc.ped"
        File kinship_tsv = cohort_prefix + "_kinship.tsv"
    }
}

task plotRelatedness {
    input {
        File kinship_tsv
        File ped_uri
        String cohort_prefix
        String plot_relatedness_script
        String hail_docker
        Int chunk_size
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(ped_uri, "GB")
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
        set -eou pipefail
        curl ~{plot_relatedness_script} > plot_relatedness.py
        python3 plot_relatedness.py ~{kinship_tsv} ~{cohort_prefix} ~{ped_uri} ~{chunk_size} > stdout
    >>>

    output {
        File relatedness_plot = "~{cohort_prefix}_relatedness_ibd0_kinship.png"
    }
}
