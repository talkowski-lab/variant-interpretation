version 1.0
    
import "Structs.wdl"

workflow deNovoSVsScatter {

    input {
        File ped_input
        File vcf_input
        File disorder_input
        String chromosome
        File raw_proband
        File raw_parents
        File raw_depth_proband
        File raw_depth_parents
        File exclude_regions
        File coverage_files
        File sample_batches
        File batch_bincov
        File python_config
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_denovo
        RuntimeAttr? runtime_attr_vcf_to_bed
    }
    
    Array[String] coverage_files = transpose(read_tsv(batch_bincov))[1]
    scatter(coverage_file in coverage_files) {
        String coverage_index_files = basename(coverage_file) + ".tbi"
    }

    # Scatter genotyping over shards
    scatter ( shard in vcf_files ) {
        call vcfToBed{
            input:
            vcf_file=shard,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_vcf_to_bed
        }

        call getDeNovo{
            input:
                bed_input=vcfToBed.bed_output,
                ped_input=ped_input,
                vcf_input=shard,
                disorder_input=disorder_input,
                chromosome=chromosome,
                raw_proband=raw_proband,
                raw_parents=raw_parents,
                raw_depth_proband=raw_depth_proband,
                raw_depth_parents=raw_depth_proband,
                exclude_regions = exclude_regions,
                coverage_files = coverage_files,
                coverage_indeces = coverage_index_files,
                sample_batches = sample_batches,
                batch_bincov = batch_bincov,
                python_config=python_config,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_denovo
        }
    }

    output {
        Array[File] per_shard_de_novo_output = runDeNovo.denovo_output
        Array[File] per_shard_de_novo_outliers = runDeNovo.denovo_outliers
        Array[File] filtered_out = runDeNovo.filtered_out
        Array[File] size_file_out = runDeNovo.size_file_out
        Array[File] coverage_output_file = runDeNovo.coverage_output_file
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
        File denovo_output = "~{basename}.denovo.bed"
        File denovo_outliers = "~{basename}.denovo.outliers.bed"
        File filtered_out = "~{basename}.filtered.txt"
        File size_file_out = "~{basename}.size.txt"
        File coverage_output_file = "~{basename}.coverage.txt"
    }

    String basename = basename(vcf_input, ".vcf.gz")
    command <<<

            bcftools view ~{vcf_input} | grep -v ^## | bgzip -c > ~{basename}.noheader.vcf.gz
            python3.9 /src/variant-interpretation/scripts/deNovoSVs.py \
                --bed ~{bed_input} \
                --ped ~{ped_input} \
                --vcf ~{basename}.noheader.vcf.gz \
                --disorder ~{disorder_input} \
                --out ~{basename}.denovo.bed \
                --raw_proband ~{raw_proband} \
                --raw_parents ~{raw_parents} \
                --raw_depth_proband ~{raw_depth_proband} \
                --raw_depth_parents ~{raw_depth_parents} \
                --config ~{python_config} \
                --filtered ~{basename}.filtered.txt \
                --size_file ~{basename}.size.txt \
                --coverage_output_file ~{basename}.coverage.txt \
                --exclude_regions ~{exclude_regions} \
                --coverage ~{batch_bincov} \
                --sample_batches ~{sample_batches} \
                --verbose True \
                --outliers ~{basename}.denovo.outliers.bed
            
            bgzip ~{basename}.denovo.bed
            bgzip ~{basename}.denovo.outliers.bed
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
