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
    File regions

    String docker_pe_evidence

    RuntimeAttr? runtime_attr_subset_sample_roi
    RuntimeAttr? runtime_attr_subset_pe_evidence


  }

  #read file and get list of samples, then:
  scatter (sample in samples) {

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

    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
    #bash automotize_PE.sh ~{sample_bed} ~{sample_batch} ~{batches_pe} > ~{sample}.pe.bed.gz

    n=0
    for i in `awk -v OFS="\t" '{print $5}' ~{sample_bed}`
    do
      n=$((n+1))
      BP1=$(awk -v OFS="\t" -v var="$n" 'NR==var {print $1":"$2}' ~{sample_bed})
      BP1_int=$(awk -v OFS="\t" -v var="$n" 'NR==var {print $1":"$2-1000"-"$2+1000}' ~{sample_bed})
      CHR2=$(awk -v OFS="\t" -v var="$n" 'NR==var {print $3}' ~{sample_bed})
      BP2=$(awk -v OFS="\t" -v var="$n" 'NR==var {print $3":"$4}' ~{sample_bed})
        # batchid=$(zcat ~{sample_batch} | awk -v OFS="\t" -v var="$i" '{if($2 == var) print $(NF)}')
        # this version is used once we don't know the exact sampleId
        # batchid=$(zcat ~{sample_batch} | fgrep -w $i | awk -v OFS="\t" '{print $(NF)}')
      batchid=$(cat ~{sample_batch} | fgrep -w $i | awk -v OFS="\t" '{print $(NF)}')
      batchfile=$(cat ~{batches_pe} | fgrep -w $batchid | awk '{print $9}')
      tabix $batchfile $BP1_int | fgrep $CHR2 > PE_metrics
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' Sample_id $i batchid $batchid batchfile $batchfile $BP1 $BP2
      cat PE_metrics
    done > ~{sample}.pe.bed.gz

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