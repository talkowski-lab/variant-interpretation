##########################
## EXPERIMENTAL WORKFLOW
##########################

# To obtains list of likely mosaic variants that failed RF due to separation only

version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcfMosaic.wdl" as MiniTasks
import "PreRFCohort.wdl" as preRF
import "MosaicDepth.wdl" as depth_mosaic
import "MosaicPesrPart1.wdl" as mosaic_pesr_part1
import "MosaicPesrPart2.wdl" as mosaic_pesr_part2

workflow MosaicManualCheck{
  input{
    File fam_file
    Int? rare_cutoff
    File outlier
    String prefix

    File sd_blacklist
    File igl_blacklist

    Array[File] per_batch_clustered_pesr_vcf_list # preRF
    Array[File] clustered_depth_vcfs
    Array[File] coverage_files
    Array[File] coverage_file_idxs
    Array[File] median_files

    Array[File] agg_metrics
    Array[File] RF_cutoffs

    String sv_pipeline_docker
    String sv_base_mini_docker
    String sv_pipeline_rdtest_docker

    RuntimeAttr? runtime_attr_concat_depth_bed
    RuntimeAttr? runtime_attr_concat_pesr_bed
    RuntimeAttr? runtime_attr_concat_depth_plot
    RuntimeAttr? runtime_attr_concat_pesr_plot

    RuntimeAttr? runtime_attr_mosaic_potential
    RuntimeAttr? runtime_attr_mosaic_depth
  }

  call subset_pesr{
      input:
          per_batch_clustered_pesr_vcf_list = per_batch_clustered_pesr_vcf_list,
          prefix = prefix,
          sv_pipeline_docker=sv_pipeline_docker
        }

  scatter (i in range(length(subset_pesr))) {
    call mosaic_pesr_part1.Mosaic as pesr1{
      input:
        name=basename(clustered_depth_vcfs[i]),
        pesr_vcfs=read_lines(subset_pesr[i]),
        metrics=agg_metrics[i],
        cutoffs=RF_cutoffs[i],
        coverage_file=coverage_files[i],
        coverage_file_idx=coverage_file_idxs[i],
        fam_file=fam_file,
        median_file=median_files[i],
        sv_pipeline_docker=sv_pipeline_docker

    }
  }
  scatter (i in range(length(clustered_depth_vcfs))) {
    call depth_mosaic.Mosaic as depth{
      input:
        name=basename(clustered_depth_vcfs[i]),
        metrics=agg_metrics[i],
        cutoffs=RF_cutoffs[i],
        rare_cutoff=rare_cutoff,
        depth_vcf=clustered_depth_vcfs[i],
        lookup=LookupGen.depthlookup,
        coverage_file=coverage_files[i],
        coverage_file_idx=coverage_file_idxs[i],
        fam_file=fam_file,
        median_file=median_files[i],
        sd_blacklist=sd_blacklist,
        igl_blacklist=igl_blacklist,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
        runtime_attr_mosaic_depth=runtime_attr_mosaic_depth,
        runtime_attr_mosaic_potential=runtime_attr_mosaic_potential
    }
  }
  call preRF.make_cohort_VCFs as LookupGen {
    input:
      pesr_vcfs = pesr1.merged_pesr,
      depth_vcfs = clustered_depth_vcfs,
      sv_pipeline_docker=sv_pipeline_docker
  }
  scatter (i in range(length(pesr1.common_potential))) {
    call mosaic_pesr_part2.Mosaic as pesr2{
      input:
        name=basename(pesr1.common_potential[i]),
        outlier=outlier,
        rare_cutoff=rare_cutoff,
        lookup=LookupGen.pesrlookup,
        potential=pesr1.common_potential[i],
        coverage_file=coverage_files[i],
        coverage_file_idx=coverage_file_idxs[i],
        fam_file=fam_file,
        median_file=median_files[i],
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker
    }
  }

  call MiniTasks.ConcatBeds as concat_depth_bed{
    input:
      shard_bed_files = depth.rare_potential,
      prefix = "~{prefix}.depth",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_depth_bed
  }

  call MiniTasks.ConcatBeds as concat_pesr_bed{
    input:
      shard_bed_files = pesr2.potentialmosaic,
      prefix = "~{prefix}.pesr",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_pesr_bed
  }
/*
  call ConcatPlots as concat_depth_plots{
    input:
      shard_plots = depth.igvplots,
      prefix = "~{prefix}.depth",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_depth_plot
  }

  call ConcatPlots as concat_pesr_plots{
    input:
      shard_plots = pesr2.igvplots,
      prefix = "~{prefix}.pesr",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_pesr_plot
  }
*/
  output{
    File depth_bed = concat_depth_bed.merged_bed_file
    File pesr_bed = concat_pesr_bed.merged_bed_file
    #File depth_plots = concat_depth_plots.merged_plots
    #File pesr_plots = concat_pesr_plots.merged_plots
  }
}

task subset_pesr{
    input {
        File per_batch_clustered_pesr_vcf_list
        String prefix
        String sv_pipeline_docker=sv_pipeline_docker
    }
    Float input_size = size(select_all([bed, ped_file]), "GB")
    Float base_mem_gb = 3.75
    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + input_size),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail
        grep "^${prefix}" ~{per_batch_clustered_pesr_vcf_list} > ~{prefix}.subset_pesr.txt
        >>>

    output{
    Array[File] subset_pesr = "~{prefix}.subset_pesr.txt"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: sv_pipeline_docker
    }

}

/*
# Merge plots from each shard:
task ConcatPlots {
  input{
    Array[File] shard_plots
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_file="~{prefix}.plot.tar.gz"

  Float input_size = size(shard_plots, "GB")
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + compression_factor * input_size,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + compression_factor)),
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

  command <<<
    set -eu

    # note head -n1 stops reading early and sends SIGPIPE to zcat,
    # so setting pipefail here would result in early termination
    mkdir output_folder/

    # no more early stopping
    set -o pipefail

    while read SPLIT; do
      tar zxf $SPLIT
      # mv plots/* output_folder/
      # argument list too long for mv
      find plots/ -type f -exec mv -t output_folder/ {} +
    done < ~{write_lines(shard_plots)} \

    tar zcf ~{prefix}.plots.tar.gz output_folder/
  >>>

  output {
    File merged_plots = "~{prefix}.plots.tar.gz"
  }
}
*/
