version 1.0
    
import "Structs2.wdl"

workflow deNovoSVsScatter {

    input {
        File ped_input
        Array[File] vcf_files
        File disorder_input
        String chromosome
        File raw_proband
        File raw_parents
        File raw_depth_proband
        File raw_depth_parents
        File exclude_regions
        File sample_batches
        File batch_bincov_index
        File python_config
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_denovo
        RuntimeAttr? runtime_attr_vcf_to_bed
        RuntimeAttr? runtime_attr_merge_bed
    }
    
    Array[String] coverage_index_files = transpose(read_tsv(batch_bincov_index))[2]

    # Scatter genotyping over shards
    scatter ( shard in vcf_files ) {
        call vcfToBed{
            input:
            vcf_file=shard,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_vcf_to_bed
        }

        call runDeNovo{
            input:
                bed_input=vcfToBed.bed_output,
                ped_input=ped_input,
                vcf_input=shard,
                disorder_input=disorder_input,
                chromosome=chromosome,
                raw_proband=raw_proband,
                raw_parents=raw_parents,
                raw_depth_proband=raw_depth_proband,
                raw_depth_parents=raw_depth_parents,
                exclude_regions = exclude_regions,
                coverage_indeces = coverage_index_files,
                sample_batches = sample_batches,
                batch_bincov_index = batch_bincov_index,
                python_config=python_config,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_denovo
        }   
    }

    call mergeBedFiles as mergeBedFilesAnnotated{
        input:
            bed_files = runDeNovo.annotation_output,
            chromosome = chromosome,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_override = runtime_attr_merge_bed
    }

    call mergeBedFiles as mergeBedFilesFinal{
        input:
            bed_files = runDeNovo.denovo_output,
            chromosome = chromosome,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_override = runtime_attr_merge_bed
    }

    output {
#        Array[File] per_shard_de_novo_output = runDeNovo.denovo_output
#        Array[File] per_shard_annotation_output = runDeNovo.annotation_output
        File per_chromosome_annotation_output_file = mergeBedFilesAnnotated.per_chromosome_denovo_output
        File per_chromosome_final_output_file = mergeBedFilesFinal.per_chromosome_denovo_output
    }
}

task runDeNovo{
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
        Array[File] coverage_indeces
        File batch_bincov_index
        File sample_batches
        File python_config
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

#    Float input_size = size(select_all([vcf_input, ped_input, disorder_input, coverage_indeces, raw_proband, raw_parents, exclude_regions, batch_bincov_index, sample_batches]), "GB")
#    Float bed_size = size(bed_input, "GB")
    Float base_mem_gb = 3.75
    Float base_disk_gb = 8

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: base_disk_gb,
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File denovo_output = "~{basename}.denovo.bed.gz"
        File annotation_output = "~{basename}.annotation.bed.gz"
    }

    String basename = basename(vcf_input, ".vcf.gz")
    command <<<
            set -euo pipefail
            bcftools view ~{vcf_input} | grep -v ^## | bgzip -c > ~{basename}.noheader.vcf.gz
            python3.9 /src/variant-interpretation/scripts/deNovoSVs.py \
                --bed ~{bed_input} \
                --ped ~{ped_input} \
                --vcf ~{basename}.noheader.vcf.gz \
                --disorder ~{disorder_input} \
                --out ~{basename}.annotation.bed \
                --out_de_novo ~{basename}.denovo.bed \
                --raw_proband ~{raw_proband} \
                --raw_parents ~{raw_parents} \
                --raw_depth_proband ~{raw_depth_proband} \
                --raw_depth_parents ~{raw_depth_parents} \
                --config ~{python_config} \
                --exclude_regions ~{exclude_regions} \
                --coverage ~{batch_bincov_index} \
                --sample_batches ~{sample_batches} \
                --verbose True
            
            bgzip ~{basename}.denovo.bed
            bgzip ~{basename}.annotation.bed
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

#    Float input_size = size(vcf_file, "GB")
    Float base_mem_gb = 3.75
    Float base_disk_gb = 8

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: base_disk_gb,
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

task mergeBedFiles{
    input{
        Array[File] bed_files
        String chromosome
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

#    Float bed_files_size = size(bed_files, "GB")
    Float base_mem_gb = 3.75
    Float base_disk_gb = 8

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: base_disk_gb,
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File per_chromosome_denovo_output = "~{chromosome}.denovo.merged.bed.gz"
    }

    command {

        zcat ${bed_files[0]} | head -n+1 > ~{chromosome}.denovo.merged.bed
        zcat ${sep=" " bed_files} | grep -v ^chrom >> ~{chromosome}.denovo.merged.bed
        bgzip ~{chromosome}.denovo.merged.bed

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