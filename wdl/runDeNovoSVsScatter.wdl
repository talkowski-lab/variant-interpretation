version 1.0
    
import "Structs.wdl"

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
    
    Array[String] coverage_files = transpose(read_tsv(batch_bincov_index))[1]
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
                coverage_files = coverage_files,
                coverage_indeces = coverage_index_files,
                sample_batches = sample_batches,
                batch_bincov_index = batch_bincov_index,
                python_config=python_config,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_denovo
        }   
    }

    call mergeBedFiles{
        input:
            bed_files = runDeNovo.denovo_output,
            chromosome = chromosome,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_override = runtime_attr_merge_bed
    }

    output {
        Array[File] per_shard_de_novo_output = runDeNovo.denovo_output
        Array[File] filtered_out = runDeNovo.filtered_out
        Array[File] size_file_out = runDeNovo.size_file_out
        Array[File] coverage_output_file = runDeNovo.coverage_output_file
        File merged_denovo_output_file = mergeBedFiles.per_chromosome_denovo_output
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
        Array[File] coverage_files
        Array[File] coverage_indeces
        File batch_bincov_index
        File sample_batches
        File python_config
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([vcf_input, bed_input, ped_input, disorder_input, raw_proband, raw_parents, exclude_regions, coverage_files, coverage_indeces, batch_bincov_index, sample_batches]), "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: ceil(base_mem_gb),
                                      disk_gb: ceil(base_disk_gb + input_size * 2.0),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File denovo_output = "~{basename}.denovo.bed.gz"
        File filtered_out = "~{basename}.filtered.txt"
        File size_file_out = "~{basename}.size.txt"
        File coverage_output_file = "~{basename}.coverage.txt"
    }

    String basename = basename(vcf_input, ".vcf.gz")
    command <<<

            bcftools view ~{vcf_input} | grep -v ^## | bgzip -c > ~{basename}.noheader.vcf.gz
            gunzip -c ~{raw_proband} > ${chromosome}.proband.reformatted.sorted.bed
            gunzip -c ~{raw_parents} > ${chromosome}.parents.reformatted.sorted.bed
            gunzip -c ~{raw_depth_proband} > ${chromosome}.proband.depth.reformatted.sorted.bed
            gunzip -c ~{raw_depth_parents} > ${chromosome}.parents.depth.reformatted.sorted.bed
            python3.9 /src/variant-interpretation/scripts/deNovoSVs.py \
                --bed ~{bed_input} \
                --ped ~{ped_input} \
                --vcf ~{basename}.noheader.vcf.gz \
                --disorder ~{disorder_input} \
                --out ~{basename}.denovo.bed \
                --raw_proband ${chromosome}.proband.reformatted.sorted.bed \
                --raw_parents ${chromosome}.parents.reformatted.sorted.bed \
                --raw_depth_proband ${chromosome}.proband.depth.reformatted.sorted.bed \
                --raw_depth_parents ${chromosome}.parents.depth.reformatted.sorted.bed \
                --config ~{python_config} \
                --filtered ~{basename}.filtered.txt \
                --size_file ~{basename}.size.txt \
                --coverage_output_file ~{basename}.coverage.txt \
                --exclude_regions ~{exclude_regions} \
                --coverage ~{batch_bincov_index} \
                --sample_batches ~{sample_batches} \
                --verbose True
            
            bgzip ~{basename}.denovo.bed
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
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(base_disk_gb + input_size * 2.0),
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

    Float bed_files_size = size(bed_files, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(base_disk_gb + (bed_files_size) * 2.0),
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
        set -eu

        zcat ${bed_files[1]} | head -n+1 > ~{chromosome}.denovo.merged.bed
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
