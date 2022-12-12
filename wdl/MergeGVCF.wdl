version 1.0

## Copyright Broad Institute, 2019
##
## This workflow merges an array of gVCFs into one
## Adapted from github.com/gatk-workflows/gatk4-germline-snps-indels/haplotypecaller-gvcf-gatk4:2.3.1
## Author: asanchis@broadinstitute.org
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.

# WORKFLOW DEFINITION
workflow HaplotypeCallerGvcf_GATK4 {
  input {
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_filename
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.0.0"
    String gatk_path = "/gatk/gatk"
  }

  # Merge per-interval GVCFs
  call MergeGVCFs {
    input:
      input_vcfs = input_vcfs,
      input_vcfs_indexes = input_vcfs_indexes,
      output_filename = output_filename,
      docker = gatk_docker,
      gatk_path = gatk_path
  }

  # Outputs that will be retained when execution is complete
  output {
    File output_vcf = MergeGVCFs.output_vcf
    File output_vcf_index = MergeGVCFs.output_vcf_index
  }
}

# TASK DEFINITIONS

# Merge GVCFs generated per-interval for the same sampls/s
task MergeGVCFs {
  input {
    # Command parameters

    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_filename

    String gatk_path

    # Runtime parameters
    String docker
    Int? mem_gb
    Int? disk_space_gb
    Int? preemptible_attempts
  }
    Boolean use_ssd = false
    Int machine_mem_gb = select_first([mem_gb, 3])
    Int command_mem_gb = machine_mem_gb - 1

  command {
  set -e

    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G"  \
      MergeVcfs \
      --INPUT ~{sep=' --INPUT ' input_vcfs} \
      --OUTPUT ~{output_filename}
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + select_first([disk_space_gb, 100]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
  }
  output {
    File output_vcf = "~{output_filename}"
    File output_vcf_index = "~{output_filename}.tbi"
  }
}

