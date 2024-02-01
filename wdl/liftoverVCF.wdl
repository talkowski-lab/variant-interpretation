version 1.0

## This workflow performs lifover of a VCF to a new genome assembly
## Author: asanchis@broadinstitute.org

# IMPORT STRUCTS
struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

# WORKFLOW DEFINITION
workflow liftoverVCF {
  input {
    File input_vcf
    File input_vcf_index
    File chain_file
    File new_reference_fasta
    File new_reference_dict

    Int java_mem

    Array[String] contigs

    String output_name
    String docker_bcftools
    String docker_gatk

    RuntimeAttr? runtime_attr_splitVCF
    RuntimeAttr? runtime_attr_liftOver
    RuntimeAttr? runtime_attr_mergeVCF

  }

  scatter (contig in contigs) {

    call splitVCF{
      input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        contig = contig,
        docker_path = docker_bcftools,
        runtime_attr_override = runtime_attr_splitVCF
      }

    call liftover {
      input:
        input_vcf = splitVCF.contig_vcf,
        input_vcf_index = splitVCF.contig_vcf_index,
        chain_file = chain_file,
        new_reference_fasta = new_reference_fasta,
        new_reference_dict = new_reference_dict,
        java_mem = java_mem,
        contig = contig,
        docker_path = docker_gatk,
        runtime_attr_override = runtime_attr_liftOver
      }
  }

  call mergeVCF{
    input:
      shard_vcf_files = liftover.contig_vcf,
      shard_vcf_files_rejected = liftover.rejected_file,
      docker_path = docker_bcftools,
      runtime_attr_override = runtime_attr_mergeVCF
    }

  output {
    File output_vcf = mergeVCF.output_vcf
    File output_vcf_index = mergeVCF.output_vcf_index
    File output_rejected = mergeVCF.output_rejected
    File output_rejected_index = mergeVCF.output_rejected_index
    }
}

# TASK DEFINITIONS
task splitVCF {
  input {
    File input_vcf
    File input_vcf_index
    String contig
    String docker_path
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output{
    File contig_vcf = "~{contig}.vcf.gz"
    File contig_vcf_index = "~{contig}.vcf.gz.tbi"
  }

  command {
    set -e
    bcftools view -r ~{contig} ~{input_vcf} -O z -o ~{contig}.vcf.gz
    tabix -p vcf ~{contig}.vcf.gz
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_path
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


task liftover {
  input {
    File input_vcf
    File input_vcf_index
    File chain_file
    File new_reference_fasta
    File new_reference_dict
    Int java_mem
    String contig
    String docker_path
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File contig_vcf = "~{contig}.lov.vcf.gz"
    File rejected_file = "~{contig}.rejected.vcf.gz"
  }

  command {
    set -euo pipefail

    gatk --java-options "-Xmx~{java_mem}g" LiftoverVcf \
      -I ~{input_vcf} \
      -O ~{contig}.lov.vcf.gz \
      --CHAIN ~{chain_file} \
      --REJECT ~{contig}.rejected.vcf.gz \
      -R ~{new_reference_fasta}
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_path
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task mergeVCF {
  input {
    Array[File] shard_vcf_files
    Array[File] shard_vcf_files_rejected
    String docker_path
    RuntimeAttr? runtime_attr_override
  }

    Float input_size = size(shard_vcf_files, "GB") + size(shard_vcf_files_rejected, "GB") 
    Float base_disk_gb = 10.0
    Float input_disk_scale = 10.0
    RuntimeAttr runtime_default = object {
        mem_gb: 8,
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
        docker: docker_path
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

  output {
    File output_vcf = "sorted_merged_lov.vcf.gz"
    File output_vcf_index = "sorted_merged_lov.vcf.gz.tbi"
    File output_rejected = "sorted_rejected.vcf.gz"
    File output_rejected_index = "sorted_rejected.vcf.gz"
  }

  command <<<
    set -euo pipefail
    # shard_vcfs=~{write_lines(shard_vcf_files)}
    # shard_vcfs_rejected=~{write_lines(shard_vcf_files_rejected)}

    # for vcf in $(cat $shard_vcfs); 
    # do
    #   tabix $vcf
    # done

    # for vcf in $(cat $shard_vcfs_rejected); 
    # do
    #   tabix $vcf
    # done

    bcftools concat -n ~{sep=" " shard_vcf_files} -O z -o merged_lov.vcf.gz
    bcftools concat -n ~{sep=" " shard_vcf_files_rejected} -O z -o rejected.vcf.gz
    
    bcftools sort -Oz -o sorted_merged_lov.vcf.gz merged_lov.vcf.gz
    bcftools sort -Oz -o sorted_rejected.vcf.gz rejected.vcf.gz

    tabix -p vcf sorted_merged_lov.vcf.gz
    tabix -p vcf sorted_rejected.vcf.gz
  >>>
}