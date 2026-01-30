version 1.0

## This workflow obtains PE or SR evidence for a list of regions in a given batch

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
workflow pe_sr_evidence {
  input {
    Array[File] batches_pe_sr
    Array[String] batch_name
    File regions
    String docker_pe_sr_evidence
    RuntimeAttr? runtime_attr_subset_pe_sr_evidence
  }

    scatter (i in range(length(batches_pe_sr))) {
      call subset_pe_sr_evidence {
        input:
          roi_bed = regions,
          batch_file = batches_pe_sr[i],
          batch_name = batch_name[i],
          docker_path = docker_pe_sr_evidence,
          runtime_attr_override = runtime_attr_subset_pe_sr_evidence
  }
}

  output {
    Array[File] final_pe_sr_evidence = subset_pe_sr_evidence.batch_pe_sr_evidence
  }

}

# TASK DEFINITIONS

task subset_pe_sr_evidence {
  input {
    File roi_bed
    File batch_file
    String batch_name
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
    File batch_pe_sr_evidence = "~{batch_name}.pe_sr.bed.gz"
  }

  command <<<
    set -ex

    grep -v '^chrom' ~{roi_bed} > roi_noheader.bed
    
    touch ~{batch_name}.pe_sr.bed

    if [[ $(wc -l <roi_noheader.bed) -ge 1 ]]; then

      export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

      while read chrom start end; do

        coords="$chrom:$start-$end"
        
        echo "Getting PE/SR evidence for ROI Coordinates: $coords in batch $batch_file"
        tabix ~{batch_file} "$coords" >> ~{batch_name}.pe_sr.bed
        
      done < roi_noheader.bed
    
    bgzip ~{batch_name}.pe_sr.bed

    fi

  >>>

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