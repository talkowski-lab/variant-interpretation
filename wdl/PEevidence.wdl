version 1.0

## This workflow obtains PE evidence for a list of samples and regions

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
workflow PEevidence {
  input {
    File batches_pe
    File sample_batch
    File sample_list
    File regions

    String docker_pe_evidence

    RuntimeAttr? runtime_attr_subset_sample_roi
    RuntimeAttr? runtime_attr_subset_pe_evidence

  }

    Array[String] samples = transpose(read_tsv(sample_list))

    scatter (sample in samples){

    call subset_sample_roi{
      input:
        sample = sample,
        regions = regions,
        docker_path = docker_pe_evidence,
        runtime_attr_override = runtime_attr_subset_sample_roi
      }

    call subset_pe_evidence {
      input:
        sample_bed = subset_sample_roi.sample_roi,
        sample_batch = sample_batch,
        batches_pe = batches_pe,
        sample = sample,
        docker_path = docker_pe_evidence,
        runtime_attr_override = runtime_attr_subset_pe_evidence
      }
  }

  output {
    Array [File] subset_sample_pe = subset_pe_evidence.sample_pe
    }

}

# TASK DEFINITIONS
task subset_sample_roi {
  input {
    File regions
    String sample
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
    File sample_roi = "~{sample}.bed.gz"
  }

  command {
    set -e
    grep -w ~{sample} ~{regions} > ~{sample}.bed.gz
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

task subset_pe_evidence {
  input {
    File sample_bed
    File sample_batch
    File batches_pe
    String sample
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
    File sample_pe = "~{sample}.pe.bed.gz"
  }

  command {
    set -euo pipefail

    i=0
    while read chr1 pos1 chr2 pos2 sample carriers svname; do
      ((i++))
      if [ "$i" -gt 1 ]; then
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
          batchname=$(grep -w $sample ~{sample_batch} | cut -f1)
          batchfile=$(grep -w $batchname ~{batches_pe} | cut -f2)
          pos1_win=$(($pos1-1000))
          pos2_win=$(($pos1+1000))
          coords="chr1:$pos1_win-$pos2_win"
          tabix $batchfile $coords | grep -w $chr2 | bgzip -c > "~{sample}.pe.bed.gz"
      else
        print("Regions for sample $sample missing")
      fi
    done < ~{sample_bed}

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