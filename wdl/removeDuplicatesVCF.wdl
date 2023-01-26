version 1.0

## This workflow removes duplicates in a VCF, required for joint genotyping
## Author: asanchis@broadinstitute.org

# WORKFLOW DEFINITION
workflow removeDuplicatesWorkflow {
  input {
    File input_vcf
    String sample_id
    String docker = "us.gcr.io/broad-gatk/gatk:4.2.0.0"
  }

  call removeDuplicatesVCF {
    input:
        File input_vcf
        String sample_id
        String docker = docker
  }

  output {
    File output_vcf = uniqueVCF.output_vcf
    File output_vcf_index = uniqueVCF.output_vcf_index
  }
}

# TASK DEFINITIONS

task removeDuplicatesVCF {
  input {

    File input_vcf
    String sample_id
    String docker

    # Runtime parameters
    Int? mem_gb
    Int? disk_space_gb
    Int? preemptible_attempts
  }
    Boolean use_ssd = false
    Int machine_mem_gb = select_first([mem_gb, 3])
    Int command_mem_gb = machine_mem_gb - 1

  command {
    set -e

    bcftools view -h ~{input_vcf} | uniq | \
        bcftools view -h ~{input_vcf} -O z -o ~{sample_id}.readgroupadded.g.vcf.gz

    bcftools index ~{sample_id}.readgroupadded.g.vcf.gz

  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + select_first([disk_space_gb, 100]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
  }
  output {
    File output_vcf = "~{sample_id}.readgroupadded.g.vcf.gz"
    File output_vcf_index = "~{sample_id}.readgroupadded.g.vcf.gz.tbi"
  }
}