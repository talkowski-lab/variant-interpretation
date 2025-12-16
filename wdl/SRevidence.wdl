version 1.0

## This workflow obtains SR evidence for a list of samples and regions

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
workflow SRevidence {
  input {
    File merged_SR
    String batch_id
    File regions
    String docker_sr_evidence
    RuntimeAttr? runtime_attr_subset_sr_evidence

  }

  call subset_sr_evidence {
    input:
      regions = regions,
      merged_SR = merged_SR,
      batch_id = batch_id,
      docker_path = docker_sr_evidence,
      runtime_attr_override = runtime_attr_subset_sr_evidence
    }

  output {
    File subset_batch_pe = subset_sr_evidence.batch_pe
    }

}

# TASK DEFINITIONS
task subset_sr_evidence {
  input {
    String batch_id
    File merged_SR
    File regions
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
    File batch_pe = "~{batch_id}.merged_SR.ROI.bed.gz"
  }

  command {
    set -ex

    export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
    tabix ~{merged_SR} ~{regions} | bgzip -c > ~{batch_id}.merged_SR.ROI.bed.gz

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