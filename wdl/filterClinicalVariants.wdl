version 1.0

import "mergeVCFs.wdl" as mergeVCFs
import "wes-denovo-helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? gpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow filterClinicalVariants {
    input {
        Array[File] annot_vcf_files
        File ped_uri
        File spliceAI_gene_annotation_file
        File pangolin_gene_annotation_file
        File ref_fasta

        String cohort_prefix
        String filter_clinical_variants_script
        String filter_clinical_variants_omim_script
        # String filter_comphets_xlr_script

        String hail_docker
        String sv_base_mini_docker
        String spliceAI_docker
        String pangolin_docker

        Float af_threshold=0.1
        Float gnomad_af_threshold=0.05
        Float am_threshold=0.56
        Float mpc_threshold=2
        Float gnomad_rec_threshold=1
        Float gnomad_dom_threshold=0.001
        Float loeuf_v2_threshold=0.35
        Float loeuf_v4_threshold=0.6

        String genome_build='GRCh38'
        Int families_per_chunk=500
        Int max_distance=50
        Boolean mask=false

        Boolean pass_filter=true
        Boolean run_pangolin=false

        RuntimeAttr? runtime_attr_merge_clinvar
        RuntimeAttr? runtime_attr_merge_omim_rec_vcfs
        RuntimeAttr? runtime_attr_merge_clinvar_vcfs
        RuntimeAttr? runtime_attr_merge_omim_dom
        RuntimeAttr? runtime_attr_merge_omim_rec
    }

    scatter (vcf_file in annot_vcf_files) {
        call runClinicalFiltering {
            input:
            vcf_file=vcf_file,
            ped_uri=ped_uri,
            filter_clinical_variants_script=filter_clinical_variants_script,
            hail_docker=hail_docker,
            af_threshold=af_threshold,
            gnomad_af_threshold=gnomad_af_threshold,
            genome_build=genome_build,
            pass_filter=pass_filter
        }

        call runClinicalFilteringOMIM {
            input:
            vcf_file=runClinicalFiltering.filtered_vcf,
            ped_uri=ped_uri,
            filter_clinical_variants_omim_script=filter_clinical_variants_omim_script,
            hail_docker=hail_docker,
            am_threshold=am_threshold,
            mpc_threshold=mpc_threshold,
            gnomad_rec_threshold=gnomad_rec_threshold,
            gnomad_dom_threshold=gnomad_dom_threshold,
            loeuf_v2_threshold=loeuf_v2_threshold,
            loeuf_v4_threshold=loeuf_v4_threshold,
            genome_build=genome_build
        }

        ## TODO: logic for spliceAI/Pangolin
        
        # if (run_pangolin) {
        #     call runPangolin {
        #         input:
        #         vcf_file=int_splice_vcf,
        #         ref_fasta=ref_fasta,
        #         gene_annotation_file=pangolin_gene_annotation_file,
        #         pangolin_docker=pangolin_docker,
        #         max_distance=max_distance,
        #         mask=mask
        #     }
        # }

        # File final_splice_vcf = select_first([runPangolin.pangolin_vcf, int_splice_vcf])
        # String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
        # String prefix = basename(vcf_file, file_ext) + '_filtered_splice_noncoding_high_moderate_impact_variants'

        # call mergeVCFs.mergeVCFs as mergeSpliceWithNonCodingImpactVars {
        #     input:  
        #         vcf_files=[final_splice_vcf, filterRareSpliceImpactVariants.noncoding_impact_vcf],
        #         sv_base_mini_docker=sv_base_mini_docker,
        #         cohort_prefix=prefix,
        #         sort_after_merge=true,
        #         naive=false   
        # }
    }   

    call helpers.mergeResultsPython as mergeClinVar {
        input:
            tsvs=runClinicalFiltering.clinvar,
            hail_docker=hail_docker,
            input_size=size(runClinicalFiltering.clinvar, 'GB'),
            merged_filename=cohort_prefix+'_clinvar_variants.tsv.gz',
            runtime_attr_override=runtime_attr_merge_clinvar
    }

    call helpers.mergeResultsPython as mergeOMIMDominant {
        input:
            tsvs=runClinicalFilteringOMIM.omim_dominant,
            hail_docker=hail_docker,
            input_size=size(runClinicalFilteringOMIM.omim_dominant, 'GB'),
            merged_filename=cohort_prefix+'_OMIM_dominant.tsv.gz',
            runtime_attr_override=runtime_attr_merge_omim_dom
    }

    call mergeVCFs.mergeVCFs as mergeOMIMRecessive {
        input:  
            vcf_files=runClinicalFilteringOMIM.omim_recessive_vcf,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix + '_OMIM_recessive',
            runtime_attr_override=runtime_attr_merge_omim_rec_vcfs
    }

    call mergeVCFs.mergeVCFs as mergeClinVarVCFs {
        input:  
            vcf_files=runClinicalFiltering.clinvar_vcf,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix + '_ClinVar_variants',
            runtime_attr_override=runtime_attr_merge_clinvar_vcfs
    }

    output {
        File clinvar_tsv = mergeClinVar.merged_tsv
        File clinvar_vcf = mergeClinVarVCFs.merged_vcf_file
        File clinvar_vcf_idx = mergeClinVarVCFs.merged_vcf_idx
        File omim_recessive_vcf = mergeOMIMRecessive.merged_vcf_file
        File omim_recessive_vcf_idx = mergeOMIMRecessive.merged_vcf_idx
        File omim_dominant_tsv = mergeOMIMDominant.merged_tsv
    }
}

task runClinicalFiltering {
    input {
        File vcf_file
        File ped_uri

        String filter_clinical_variants_script
        String hail_docker
        String genome_build

        Float af_threshold
        Float gnomad_af_threshold
        Boolean pass_filter

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
        curl ~{filter_clinical_variants_script} > filter_vcf.py
        python3 filter_vcf.py ~{vcf_file} ~{prefix} ~{cpu_cores} ~{memory} \
            ~{ped_uri} ~{af_threshold} ~{gnomad_af_threshold} ~{genome_build} ~{pass_filter}
    }

    output {
        File clinvar = prefix + '_clinvar_variants.tsv.gz'
        File clinvar_vcf = prefix + '_clinvar_variants.vcf.bgz'
        File filtered_vcf = prefix + '_clinical.vcf.bgz'
    }
}

task filterCompHetsXLR {
    input {
        File vcf_file
        File ped_uri

        String filter_comphets_xlr_script
        String hail_docker
        String genome_build 

        Float af_threshold
        Float gnomad_af_threshold
        Float am_threshold
        Float mpc_threshold
        Float gnomad_rec_threshold
        Float gnomad_dom_threshold
        Float loeuf_v2_threshold
        Float loeuf_v4_threshold

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
        curl ~{filter_comphets_xlr_script} > filter_vcf.py
        python3 filter_vcf.py ~{vcf_file} ~{prefix} ~{cpu_cores} ~{memory} \
            ~{ped_uri} ~{af_threshold} ~{gnomad_af_threshold} ~{am_threshold} \
            ~{mpc_threshold} ~{gnomad_rec_threshold} ~{gnomad_dom_threshold} \
            ~{loeuf_v2_threshold} ~{loeuf_v4_threshold} ~{genome_build}
    }

    output {
        File omim_recessive_comphet_xlr = prefix + '_OMIM_recessive_comphet_XLR.tsv.gz'
    }
}

task runClinicalFilteringOMIM {
    input {
        File vcf_file
        File ped_uri

        String filter_clinical_variants_omim_script
        String hail_docker
        String genome_build
        
        Float am_threshold
        Float mpc_threshold
        Float gnomad_rec_threshold
        Float gnomad_dom_threshold
        Float loeuf_v2_threshold
        Float loeuf_v4_threshold

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
        curl ~{filter_clinical_variants_omim_script} > filter_vcf.py
        python3 filter_vcf.py ~{vcf_file} ~{prefix} ~{cpu_cores} ~{memory} ~{ped_uri} \
            ~{am_threshold} ~{mpc_threshold} ~{gnomad_rec_threshold} ~{gnomad_dom_threshold} \
            ~{loeuf_v2_threshold} ~{loeuf_v4_threshold} ~{genome_build}
    }

    output {
        File omim_recessive_vcf = prefix + '_OMIM_recessive.vcf.bgz'
        File omim_dominant = prefix + '_OMIM_dominant.tsv.gz'
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
        tabix ~{vcf_file}
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