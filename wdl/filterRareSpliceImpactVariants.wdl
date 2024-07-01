version 1.0

import "mergeVCFs.wdl" as mergeVCFs

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? gpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow getRareSpliceImpactVariants {
    input {
        Array[File] vep_vcf_files
        File gene_annotation_file
        File ref_fasta

        String cohort_prefix
        String filter_rare_splice_impact_script
        String hail_docker
        String sv_base_mini_docker
        String spliceAI_docker
        String pangolin_docker

        Int ac_threshold=20
        Float gnomad_af_threshold=0.05
        Int max_distance=50
        Boolean mask=false

        Boolean run_spliceAI=true
        Boolean run_pangolin=false
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

        if (run_spliceAI) {
            call runSpliceAI {
                input:
                vcf_file=filterRareSpliceImpactVariants.splice_vcf,
                vcf_idx=filterRareSpliceImpactVariants.splice_vcf_idx,
                ref_fasta=ref_fasta,
                gene_annotation_file=gene_annotation_file,
                spliceAI_docker=spliceAI_docker,
                max_distance=max_distance,
                mask=mask
            }
        }
        File int_splice_vcf = select_first([runSpliceAI.spliceAI_vcf, filterRareSpliceImpactVariants.splice_vcf])
        
        if (run_pangolin) {
            call runPangolin {
                input:
                vcf_file=int_splice_vcf,
                vcf_idx=int_splice_vcf+'.tbi',
                ref_fasta=ref_fasta,
                gene_annotation_file=gene_annotation_file,
                pangolin_docker=pangolin_docker,
                max_distance=max_distance,
                mask=mask
            }
        }

        File final_splice_vcf = select_first([runPangolin.pangolin_vcf, int_splice_vcf])
        String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
        String prefix = basename(vcf_file, file_ext) + '_filtered_splice_noncoding_high_moderate_impact_variants'

        call mergeVCFs.mergeVCFs as mergeSpliceWithNonCodingImpactVars {
            input:  
                vcf_files=[final_splice_vcf, filterRareSpliceImpactVariants.noncoding_impact_vcf],
                sv_base_mini_docker=sv_base_mini_docker,
                cohort_prefix=cohort_prefix + '_filtered_splice_noncoding_high_moderate_impact_variants',
                sort_after_merge=true        
        }
    }   

    output {
        Array[File] clinvar_vcf_files = filterRareSpliceImpactVariants.clinvar_vcf
        Array[File] clinvar_vcf_idx = filterRareSpliceImpactVariants.clinvar_vcf_idx
        Array[File] splice_nc_impact_vcf_files = mergeSpliceWithNonCodingImpactVars.merged_vcf_file
        Array[File] splice_nc_impact_vcf_idx = mergeSpliceWithNonCodingImpactVars.merged_vcf_idx
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
    String prefix = basename(vcf_file, file_ext) + '_filtered'

    command {
        curl ~{filter_rare_splice_impact_script} > filter_vcf.py
        python3 filter_vcf.py ~{vcf_file} ~{prefix} ~{cpu_cores} ~{memory} ~{ac_threshold} ~{gnomad_af_threshold}
    }

    output {
        File clinvar_vcf = prefix + '_clinvar_variants.vcf.bgz'
        File clinvar_vcf_idx = clinvar_vcf + '.tbi'
        File splice_vcf = prefix + '_splice_variants.vcf.bgz'
        File splice_vcf_idx = splice_vcf + '.tbi'
        File noncoding_impact_vcf = prefix + '_noncoding_high_moderate_impact_variants.vcf.bgz'
        File noncoding_impact_vcf_idx = noncoding_impact_vcf + '.tbi'
    }
}

task getNonEmptyVCFs {
    input {
        Array[File] vcf_files
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_files, 'GB')
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
    set -eou pipefail
    mkdir output_vcfs
    for vcf_file in $(cat ~{write_lines(vcf_files)});
        do
        tabix $vcf_file
        if [[ $(bcftools index -n $vcf_file".tbi") != 0 ]]; then
            mv $vcf_file output_vcfs/
            mv $vcf_file".tbi" output_vcfs/
        fi
    done    
    >>>

    output {
        Array[File] output_vcfs = glob('output_vcfs/*.vcf.bgz')
        Array[File] output_vcfs_idx = glob('output_vcfs/*.vcf.bgz.tbi')
    }
}

task runSpliceAI {
    input {
        File vcf_file
        File vcf_idx
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
        gpu_cores: 2,
        cpu_cores: 1, 
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    Int gpu_cores = select_first([runtime_override.gpu_cores, runtime_default.gpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        gpuType: "nvidia-tesla-t4"
        gpuCount: gpu_cores
        nvidiaDriverVersion: "450.80.02"
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: spliceAI_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String mask_str = if mask then '--mask' else ''
    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    String output_filename = basename(vcf_file, file_ext) + '_spliceAI.vcf'

    command {
        set -eou pipefail
        spliceai.py -r ~{ref_fasta} -a ~{gene_annotation_file} -i ~{vcf_file} -o ~{output_filename} \
            -d ~{max_distance} ~{mask_str} --preprocessing_threads 4
        bgzip ~{output_filename}
    }

    output {
        File spliceAI_vcf = output_filename + '.gz'
    }
}


task runPangolin {
    input {
        File vcf_file
        File vcf_idx
        File gene_annotation_file
        File ref_fasta

        Boolean mask
        Int max_distance  # Maximum distance between the variant and gained/lost splice site. An integer in the range [0, 5000]. Defaults to 50.
        String pangolin_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_file, 'GB')
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        gpu_cores: 2,
        cpu_cores: 1, 
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    Int gpu_cores = select_first([runtime_override.gpu_cores, runtime_default.gpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        gpuType: "nvidia-tesla-t4"
        gpuCount: gpu_cores
        nvidiaDriverVersion: "450.80.02"
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: pangolin_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String mask_str = if mask then '-m True' else ''
    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    String output_filename = basename(vcf_file, file_ext) + '_pangolin'

    command {
        set -eou pipefail
        zcat ~{vcf_file} > ~{basename(vcf_file, file_ext)}.vcf
        pangolin ~{basename(vcf_file, file_ext)}.vcf ~{ref_fasta} ~{gene_annotation_file} ~{output_filename} \
            -d ~{max_distance} ~{mask_str} 
        bgzip ~{output_filename}.vcf
    }

    output {
        File pangolin_vcf = output_filename + '.vcf.gz'
    }
}