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
        File raw_files_list
        File depth_raw_files_list
        File exclude_regions
        File sample_batches
        File batch_bincov
        Int records_per_shard
        String prefix
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
    
    }

    call raw.reformatRawFiles as reformatRawFiles {
        input:
            contigs = contigs,
            raw_files_list = raw_files_list,
            ped_input = ped_input,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_vcf_to_bed = runtime_attr_raw_vcf_to_bed,
            runtime_attr_merge_bed = runtime_attr_raw_merge_bed,
            runtime_attr_divide_by_chrom = runtime_attr_raw_divide_by_chrom,
            runtime_attr_reformat_bed = runtime_attr_raw_reformat_bed
    }

    call raw.reformatRawFiles as reformatDepthRawFiles {
        input:
            contigs = contigs,
            raw_files_list = depth_raw_files_list,
            ped_input = ped_input,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_vcf_to_bed = runtime_attr_raw_vcf_to_bed,
            runtime_attr_merge_bed = runtime_attr_raw_merge_bed,
            runtime_attr_divide_by_chrom = runtime_attr_raw_divide_by_chrom,
            runtime_attr_reformat_bed = runtime_attr_raw_reformat_bed
    }
      
    call getGenomicDisorders{
        input:
            genomic_disorder_input=genomic_disorder_input,
            vcf_file=vcf_file,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_gd
    }
    
    scatter (i in range(length(contigs))){
        call subsetVcf {
            input:
                vcf_file=vcf_file,
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
                ped_input=ped_input,
                vcf_input=SplitVcf.shards,
                disorder_input=getGenomicDisorders.gd_output,
                chromosome=contigs[i],
                raw_proband=reformatRawFiles.reformatted_proband_raw_files[i],
                raw_parents=reformatRawFiles.reformatted_parents_raw_files[i],
                raw_depth_proband=reformatDepthRawFiles.reformatted_proband_raw_files[i],
                raw_depth_parents=reformatDepthRawFiles.reformatted_parents_raw_files[i],
                exclude_regions = exclude_regions,
                sample_batches = sample_batches,
                batch_bincov = batch_bincov,
                python_config=python_config,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_denovo = runtime_attr_denovo,
                runtime_attr_vcf_to_bed = runtime_attr_vcf_to_bed
        }
    }

    call plot_mergeFinalBedFiles{
        input:
            bed_files = flatten(getDeNovo.per_shard_de_novo_output),
            outliers_files = flatten(getDeNovo.per_shard_de_novo_outliers),
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_merge_final_bed_files
    }

    call plot_createPlots{
        input:
            bed_file = plot_mergeFinalBedFiles.final_denovo_output,
            outliers_file = plot_mergeFinalBedFiles.final_denovo_outliers_output,
            ped_input = ped_input,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_create_plots
    }

    output {
    
        File denovo_output = plot_mergeFinalBedFiles.final_denovo_output
        File denovo_outliers_output = plot_mergeFinalBedFiles.final_denovo_outliers_output
        File denovo_output_plots = plot_createPlots.output_plots
        Array[Array[File]] filtered_out = getDeNovo.filtered_out
        Array[Array[File]] size_file_out = getDeNovo.size_file_out
        Array[Array[File]] coverage_file_out = getDeNovo.coverage_output_file
    }
}

task getDeNovo{
    input{
        File bed_input
        File ped_input
        File vcf_input
        File disorder_input
        String chromosome
        File raw_proband
        File raw_parents
        File raw_depth_proband
        File raw_depth_parents
        File exclude_regions
        Array[File] coverage_files
        Array[File] coverage_indeces
        File batch_bincov
        File sample_batches
        File python_config
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([vcf_input, bed_input, ped_input, disorder_input, raw_proband, raw_parents, exclude_regions, coverage_files, coverage_indeces, batch_bincov, sample_batches]), "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: ceil(base_mem_gb + input_size * 3.0),
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File denovo_output = "~{chromosome}.denovo.bed"
        File denovo_outliers = "~{chromosome}.denovo.outliers.bed"
        File filtered_out = "~{chromosome}.filtered.txt"
        File size_file_out = "~{chromosome}.size.txt"
    }

    String basename = basename(vcf_input, ".vcf.gz")
    command <<<

            bcftools view ~{vcf_input} | grep -v ^## | bgzip -c > ~{basename}.noheader.vcf.gz
            python3.9 /src/variant-interpretation/scripts/deNovoSVs.py \
                --bed ~{bed_input} \
                --ped ~{ped_input} \
                --vcf ~{basename}.noheader.vcf.gz \
                --disorder ~{disorder_input} \
                --out ~{chromosome}.denovo.bed \
                --raw_proband ~{raw_proband} \
                --raw_parents ~{raw_parents} \
                --raw_depth_proband ~{raw_depth_proband} \
                --raw_depth_parents ~{raw_depth_parents} \
                --config ~{python_config} \
                --filtered ~{chromosome}.filtered.txt \
                --size_file ~{chromosome}.size.txt \
                --coverage_output_file ~{chromosome}.coverage.txt \
                --exclude_regions ~{exclude_regions} \
                --coverage ~{batch_bincov} \
                --sample_batches ~{sample_batches} \
                --verbose True \
                --outliers ~{chromosome}.denovo.outliers.bed
            
            bgzip ~{chromosome}.denovo.bed
            bgzip ~{chromosome}.denovo.outliers.bed
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

task vcfToBed{
    input{
        File vcf_file
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
        File bed_output = "~{basename}.bed.gz"
    }

    String basename = basename(vcf_file, ".vcf.gz")
    command <<<
        set -euo pipefail

        svtk vcf2bed ~{vcf_file} --info ALL --include-filters ~{basename}.bed
        bgzip ~{basename}.bed
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
        Array[File] outliers_files
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float bed_files_size = size(bed_files, "GB")
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
        File final_denovo_outliers_output = "final.denovo.outliers.merged.bed.gz"
    }

    command {
        set -euo pipefail

        zcat ${bed_files[1]} | head -n+1 > final.denovo.merged.bed
        zcat ${sep=" " bed_files} | tail -n+2 -q >> final.denovo.merged.bed
        bgzip final.denovo.merged.bed
        
        head -n+1 ${outliers_files[1]} > final.denovo.outliers.merged.bed
        tail -n+2 -q ${sep=" " outliers_files} >> final.denovo.outliers.merged.bed
        bgzip final.denovo.outliers.merged.bed
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
