version 1.0

import "mergeVCFs.wdl" as mergeVCFs

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow filterRunSpliceAI {
    input {
        Array[File] vep_vcf_files
        File gene_annotation_file
        File ref_fasta

        String cohort_prefix
        String filter_rare_splice_impact_script
        String hail_docker
        String sv_base_mini_docker
        String spliceAI_docker

        Int ac_threshold=20
        Float gnomad_af_threshold=0.05
        Int max_distance=50
        Boolean mask=false
        Boolean sort_after_merge=false
    }
    
    scatter (vcf_file in vep_vcf_files) {
        call filterRareSpliceImpactVariants {
            input:
            vcf_file=vcf_file,
            filter_rare_splice_impact_script=filter_rare_splice_impact_script,
            hail_docker=hail_docker,
            ac_threshold=ac_threshold,
            gnomad_af_threshold=gnomad_af_threshold
        }

        call runSpliceAI {
            input:
            vcf_file=filterRareSpliceImpactVariants.filtered_vcf,
            ref_fasta=ref_fasta,
            gene_annotation_file=gene_annotation_file,
            spliceAI_docker=spliceAI_docker,
            max_distance=max_distance,
            mask=mask
        }
    }

    call mergeVCFs.mergeVCFs as mergeVCFs {
        input:  
            vcf_files=runSpliceAI.spliceAI_vcf,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix + '_filtered_spliceAI',
            sort_after_merge=sort_after_merge        
    }

    output {
        File filtered_spliceAI_vcf = mergeVCFs.merged_vcf_file
    }
}

task filterRareSpliceImpactVariants {
    input {
        File vcf_file
        String filter_rare_splice_impact_script
        String hail_docker
        Int ac_threshold
        Float gnomad_af_threshold
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_file, 'GB')
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

    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    String output_filename = basename(vcf_file, file_ext) + '_filtered' + file_ext

    command {
        curl ~{filter_rare_splice_impact_script} > filter_vcf.py
        python3 filter_vcf.py ~{vcf_file} ~{output_filename} ~{cpu_cores} ~{memory} ~{ac_threshold} ~{gnomad_af_threshold}
    }

    output {
        File filtered_vcf = output_filename
    }
}

task runSpliceAI {
    input {
        File vcf_file
        File gene_annotation_file
        File ref_fasta

        Boolean mask
        Int max_distance  # Maximum distance between the variant and gained/lost splice site. An integer in the range [0, 5000]. Defaults to 50.
        String spliceAI_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_file, 'GB')
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
        docker: spliceAI_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
        gpuType: "nvidia-tesla-k80"
        gpuCount: 2
        nvidiaDriverVersion: "418.87.00"
    }

    String mask_str = if mask then '--mask' else ''
    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    String output_filename = basename(vcf_file, file_ext) + '_spliceAI' + file_ext

    command {
        set -eou pipefail
        spliceai.py -r ~{ref_fasta} -a ~{gene_annotation_file} -i ~{vcf_file} -o ~{output_filename} \
            -d ~{max_distance} ~{mask_str} --preprocessing_threads 4
    }

    output {
        File spliceAI_vcf = output_filename
    }
}