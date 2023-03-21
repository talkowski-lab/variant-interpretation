version 1.0
    
import "Structs.wdl"
import "reformatRawFiles.wdl" as raw
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "runDeNovoSVsScatter.wdl" as runDeNovo

workflow deNovoSV {

    input {

        File ped_input
        File python_config
        File vcf_file
        Array[String] contigs
        File genomic_disorder_input
        File batch_raw_file
        File batch_depth_raw_file
        File exclude_regions
        File sample_batches
        File batch_bincov_index
        Int records_per_shard
        String prefix
        String? fam_ids
        String variant_interpretation_docker
        String sv_pipeline_updates_docker
        RuntimeAttr? runtime_attr_gd
        RuntimeAttr? runtime_attr_denovo
        RuntimeAttr? runtime_attr_raw_vcf_to_bed
        RuntimeAttr? runtime_attr_raw_merge_bed
        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_vcf_to_bed
        RuntimeAttr? runtime_attr_raw_divide_by_chrom
        RuntimeAttr? runtime_attr_raw_reformat_bed
        RuntimeAttr? runtime_attr_merge_final_bed_files
        RuntimeAttr? runtime_attr_create_plots
        RuntimeAttr? runtime_override_shard_vcf
        RuntimeAttr? runtime_attr_clean_ped
        RuntimeAttr? runtime_attr_call_outliers
        RuntimeAttr? runtime_attr_get_raw_files
    
    }

    if (defined(fam_id)){
        String fam_ids_ = select_first([fam_ids])
        call getBatchedFiles{
            input:
                batch_raw_file = batch_raw_file,
                batch_depth_raw_file = batch_depth_raw_file,
                fam_ids = fam_ids,
                vcf_file = vcf_file,
                sample_batches = sample_batches,
                batch_bincov_index = batch_bincov_index,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_get_raw_files
        }
    }

    call cleanPed{
        input:
            ped_input = ped_input,
            vcf_input = select_first([getBatchedFiles.subset_vcf, vcf_file]),
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_clean_ped
    }

    call raw.reformatRawFiles as reformatRawFiles {
        input:
            contigs = contigs,
            raw_files_list = select_first([getBatchedFiles.batch_raw_files_list, batch_raw_file]),
            ped_input = cleanPed.cleaned_ped,
            depth = false,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_vcf_to_bed = runtime_attr_raw_vcf_to_bed,
            runtime_attr_merge_bed = runtime_attr_raw_merge_bed,
            runtime_attr_divide_by_chrom = runtime_attr_raw_divide_by_chrom,
            runtime_attr_reformat_bed = runtime_attr_raw_reformat_bed
    }

    call raw.reformatRawFiles as reformatDepthRawFiles {
        input:
            contigs = contigs,
            raw_files_list = select_first([getBatchedFiles.batch_depth_raw_files_list, batch_depth_raw_file]),
            ped_input = cleanPed.cleaned_ped,
            depth = true,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_vcf_to_bed = runtime_attr_raw_vcf_to_bed,
            runtime_attr_merge_bed = runtime_attr_raw_merge_bed,
            runtime_attr_divide_by_chrom = runtime_attr_raw_divide_by_chrom,
            runtime_attr_reformat_bed = runtime_attr_raw_reformat_bed
    }
      
    call getGenomicDisorders{
        input:
            genomic_disorder_input=genomic_disorder_input,
            vcf_file = select_first([getBatchedFiles.subset_vcf, vcf_file]),
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_gd
    }
    
    scatter (i in range(length(contigs))){
        call subsetVcf {
            input:
                vcf_file = select_first([getBatchedFiles.subset_vcf, vcf_file]),
                chromosome=contigs[i],
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call MiniTasks.ScatterVcf as SplitVcf {
            input:
                vcf=subsetVcf.vcf_output,
                prefix=prefix,
                records_per_shard=records_per_shard,
                sv_pipeline_docker=sv_pipeline_updates_docker,
                runtime_attr_override=runtime_override_shard_vcf
        }
    
        call runDeNovo.deNovoSVsScatter as getDeNovo {
            input:
                ped_input=cleanPed.cleaned_ped,
                vcf_files=SplitVcf.shards,
                disorder_input=getGenomicDisorders.gd_output,
                chromosome=contigs[i],
                raw_proband=reformatRawFiles.reformatted_proband_raw_files[i],
                raw_parents=reformatRawFiles.reformatted_parents_raw_files[i],
                raw_depth_proband=reformatDepthRawFiles.reformatted_proband_raw_files[i],
                raw_depth_parents=reformatDepthRawFiles.reformatted_parents_raw_files[i],
                exclude_regions = exclude_regions,
                sample_batches = sample_batches,
                batch_bincov_index = select_first([getBatchedFiles.batch_bincov_index_subset, batch_bincov_index]),
                python_config=python_config,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_denovo = runtime_attr_denovo,
                runtime_attr_vcf_to_bed = runtime_attr_vcf_to_bed
        }
    }

    call plot_mergeFinalBedFiles{
        input:
            bed_files = getDeNovo.merged_denovo_output_file,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_merge_final_bed_files
    }

    call callOutliers {
        input:
            bed_file = plot_mergeFinalBedFiles.merged_output,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_call_outliers
    }

    call plot_createPlots{
        input:
            bed_file = callOutliers.final_denovo_output,
            outliers_file = callOutliers.final_denovo_outliers_output,
            ped_input = ped_input,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_create_plots
    }

    output {
    
        File denovo_output = callOutliers.final_denovo_output
        File denovo_outliers_output = callOutliers.final_denovo_outliers_output
        File denovo_output_plots = plot_createPlots.output_plots
        File per_chrom_plot = plot_createPlots.per_chrom_plot
        File per_sample_plot = plot_createPlots.per_sample_plot
        File per_freq_plot = plot_createPlots.per_freq_plot
        File per_freq_gd_plot = plot_createPlots.per_freq_gd_plot
        File per_freq_not_gd = plot_createPlots.per_freq_not_gd
        File size_plot = plot_createPlots.size_plot
        File evidence_plot = plot_createPlots.evidence_plot
        File annotation_plot = plot_createPlots.annotation_plot
        File per_type_plot = plot_createPlots.per_type_plot
        Array[Array[File]] filtered_out = getDeNovo.filtered_out
        Array[Array[File]] size_file_out = getDeNovo.size_file_out
        Array[Array[File]] coverage_file_out = getDeNovo.coverage_output_file
    }
}

task subsetVcf{
    input{
        File vcf_file
        String chromosome
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb + input_size * 3.0,
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File vcf_output = "~{chromosome}.vcf.gz"
        File no_header_vcf_output = "~{chromosome}.noheader.vcf.gz"
    }

    command <<<
        set -euo pipefail

        bcftools index ~{vcf_file}

        bcftools view ~{vcf_file} --regions ~{chromosome} -O z -o  ~{chromosome}.vcf.gz

        bcftools view ~{chromosome}.vcf.gz | grep -v ^## | bgzip -c > ~{chromosome}.noheader.vcf.gz

        
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task getGenomicDisorders{
    input{
        File vcf_file
        File genomic_disorder_input
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb + input_size * 3.0,
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File gd_output = "annotated.gd.variants.names.txt"
    }

    command <<<
        set -euo pipefail

        bedtools intersect -wa -wb -f 0.3 -r -a ~{vcf_file} -b ~{genomic_disorder_input} | cut -f 3 |sort -u> annotated.gd.variants.names.txt
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task plot_mergeFinalBedFiles{
    input{
        Array[File] bed_files
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float bed_files_size = size(bed_files, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb + (bed_files_size) * 3.0,
                                      disk_gb: ceil(base_disk_gb + (bed_files_size) * 5.0),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File merged_output = "merged.bed.gz"
    }

    command {
        set -eu

        zcat ${bed_files[1]} | head -n+1 > final.denovo.merged.bed
        zcat ${sep=" " bed_files} | grep -v ^chrom >> final.denovo.merged.bed
        bgzip final.denovo.merged.bed
    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task callOutliers{
    input{
        File bed_file
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float bed_files_size = size(bed_file, "GB")
    Float outliers_files_size = size(outliers_files, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb + (bed_files_size + outliers_files_size) * 3.0,
                                      disk_gb: ceil(base_disk_gb + (bed_files_size + outliers_files_size) * 5.0),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File final_denovo_output = "final.denovo.merged.bed.gz"
        File final_denovo_outliers_output = "final.denovo.merged.outliers.bed.gz"
    }

    command {
        set -euo pipefail

        python3.9 /src/variant-interpretation/scripts/deNovoOutliers.py --bed ~{bed_file}
        bgzip final.denovo.merged.bed
        bgzip final.denovo.merged.outliers.bed

    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task plot_createPlots{
    input{
        File bed_file
        File outliers_file
        File ped_input
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([bed_file, outliers_file]), "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb + input_size * 3.0,
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File output_plots = "output_plots.pdf"
        File per_chrom_plot = "per_chrom.png"
        File per_sample_plot = "per_sample.png"
        File per_freq_plot = "per_freq.png"
        File per_freq_gd_plot = "per_freq_gd.png"
        File per_freq_not_gd = "per_freq_not_gd.png"
        File size_plot = "size.png"
        File evidence_plot = "evidence.png"
        File annotation_plot = "annotation.png"
        File per_type_plot = "per_type.png"

    }

    command {
        set -euo pipefail

        Rscript /src/variant-interpretation/scripts/denovoSV_plots.R ${bed_file} ${outliers_file} ${ped_input} output_plots.pdf

    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task cleanPed{
    input{
        File ped_input
        File vcf_input
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([ped_input, vcf_input]), "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb + input_size * 3.0,
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File cleaned_ped = "subset_cleaned_ped.txt"
    }

    command {
        set -euo pipefail

        Rscript /src/variant-interpretation/scripts/cleanPed.R ${ped_input}
        cut -f2 cleaned_ped.txt | tail -n+2 > all_samples.txt
        bcftools query -l ${vcf_input} > samples_to_include_in_ped.txt
        grep -w -v -f samples_to_include_in_ped.txt all_samples.txt > excluded_samples.txt
        grep -w -f excluded_samples.txt cleaned_ped.txt | cut -f1 | sort -u > excluded_families.txt
        grep -w -v -f excluded_families.txt cleaned_ped.txt > subset_cleaned_ped.txt

    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task getBatchedFiles{
    input{
        File batch_raw_file
        File batch_depth_raw_file
        File fam_ids
        File ped_input
        File sample_batches
        File batch_bincov_index
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([ped_input, vcf_input]), "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb + input_size * 3.0,
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File batch_raw_files_list = "batch_raw_files_list.txt"
        File batch_depth_raw_files_list = "batch_depth_raw_files_list.txt"
        File batch_bincov_index_subset = "batch_bincov_index.txt"
        File subset_vcf = "subset_vcf_file.vcf"
    }

    command {
        set -euo pipefail
        grep -w -f ${fam_id} ${ped_input} | cut -f2 | sort -u > samples.txt
        grep -w -f samples.txt ${sample_batches} | cut -f2 | sort -u > batches.txt
        grep -w -f batches.txt ${batch_bincov_index} > batch_bincov_index.txt
        grep -w -f batches.txt ${batch_raw_file} > batch_raw_files_list.txt
        grep -w -f batches.txt ${batch_depth_raw_file} > batch_depth_raw_files_list.txt
        grep -w -f samples.txt ${vcf_file} > subset_vcf_file.vcf
    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}
