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
#    File sample_batch
#    File sample_list
    File regions

    String docker_pe_evidence

    RuntimeAttr? runtime_attr_subset_tloc_roi
    RuntimeAttr? runtime_attr_subset_pe_evidence

  }

    Array[String] tloc_ids = transpose(read_tsv(regions))[4]

    scatter (tloc in tloc_ids){

      call subset_tloc_roi{
        input:
          name = tloc,
          regions = regions,
          docker_path = docker_pe_evidence,
          runtime_attr_override = runtime_attr_subset_tloc_roi
      }

      call subset_pe_evidence {
        input:
          tloc_bed = subset_tloc_roi.tloc_roi,
  #        sample_batch = sample_batch,
          batches_pe = batches_pe,
          name = tloc,
          docker_path = docker_pe_evidence,
          runtime_attr_override = runtime_attr_subset_pe_evidence
      }
    }

  output {
    Array [File] subset_pe_evidence = subset_pe_evidence.tloc_pe
  }

}

# TASK DEFINITIONS
task subset_tloc_roi {
  input {
    File regions
    String name
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
    File tloc_roi = "~{name}.bed.gz"
  }

  command {
    set -euo pipefail

    zcat ~{regions} | \
    grep -Ew "^chrom1|~{name}" | bgzip -c > ~{name}.bed.gz
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
    File tloc_bed
#    File sample_batch
    File batches_pe
    String name
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
    File tloc_pe = "~{name}.pe.bed.gz"
  }

  command <<<
    set -ex
    zcat ~{tloc_bed} | grep -v ^chrom > tloc_noheader.pe.bed


    if [[ $(wc -l <tloc_noheader.pe.bed) -ge 1 ]]; then
      while read chr1 pos1 chr2 pos2 svname; do
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        pos1_win=$(($pos1-2000))
        pos2_win=$(($pos1+2000))
        coords="$chr1:$pos1_win-$pos2_win"

        awk -v OFS="\t" '{if (NR>1) print $2}' ~{batches_pe} | while read batch_id; do
          echo "Getting PE evidence for Tloc $svname Coordinates: $coords in batch $batch_id"
          tabix $batchfile $coords | grep -w $chr2 >> ~{name}.pe.bed
        done
      done < tloc_noheader.pe.bed

      bgzip ~{name}.pe.bed

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