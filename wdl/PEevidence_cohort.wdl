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
    File regions
    File regions_sample

    String docker_pe_evidence

    RuntimeAttr? runtime_attr_subset_tloc_roi
    RuntimeAttr? runtime_attr_subset_pe_evidence
    RuntimeAttr? runtime_attr_calculate_freq
    RuntimeAttr? runtime_attr_merge_freq
    RuntimeAttr? runtime_attr_summary
    RuntimeAttr? runtime_attr_merge_summary

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

      call calculate_frequencies {
        input:
          tloc_bed = subset_tloc_roi.tloc_roi,
          tloc_pe_evidence = subset_pe_evidence.tloc_pe,
          name = tloc,
          docker_path = docker_pe_evidence,
          runtime_attr_override = runtime_attr_calculate_freq
      }

      call create_summary {
        input:
          tloc_sample_bed = regions_sample,
          tloc_info_file = calculate_frequencies.tloc_info,
          name = tloc,
          docker_path = docker_pe_evidence,
          runtime_attr_override = runtime_attr_summary
      }
    }

    call merge_frequencies {
      input:
        tloc_freq_files = calculate_frequencies.tloc_freq,
        docker_path = docker_pe_evidence,
        runtime_attr_override = runtime_attr_merge_freq
    }

    call merge_summary {
      input:
        tloc_summary_files = create_summary.tloc_summary,
        docker_path = docker_pe_evidence,
        runtime_attr_override = runtime_attr_merge_summary
    }

  output {
    Array [File] tlocs_pe_evidence = subset_pe_evidence.tloc_pe
    Array [File] tlocs_info = calculate_frequencies.tloc_info
    File tlocs_frequencies = merge_frequencies.tloc_merged_freq
    File tlocs_summary = merge_summary.tloc_merged_summary
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
    File tloc_roi = "~{name}.bed"
  }

  command {
    set -euo pipefail

    grep -Ew "^chrom1|~{name}" ~{regions} > ~{name}.bed
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
    grep -v ^chrom ~{tloc_bed} > tloc_noheader.pe.bed


    if [[ $(wc -l <tloc_noheader.pe.bed) -ge 1 ]]; then
      while read chr1 pos1 chr2 pos2 svname; do
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        pos1_win=$(($pos1-2000))
        pos2_win=$(($pos1+2000))
        coords="$chr1:$pos1_win-$pos2_win"

        awk -v OFS="\t" '{if (NR>1) print $2}' ~{batches_pe} | while read batch_id; do
          echo "Getting PE evidence for Tloc $svname Coordinates: $coords in batch $batch_id"
          tabix $batch_id $coords | grep -w $chr2 >> ~{name}.pe.bed
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

task calculate_frequencies {
  input {
    File tloc_bed
    File tloc_pe_evidence
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
    File tloc_freq ="~{name}.freq"
    File tloc_info ="~{name}.info"
  }

  command <<<
    set -ex

    i=1
    cat ~{tloc_bed} | while read line
    do
    i=$((i+1))
    BP1_interval=$(echo $line | awk -v OFS="\t" '{print $1":"$2-2000"-"$2+2000}')
    BP2_interval=$(echo $line | awk -v OFS="\t" '{print $3":"$4}')
    BP2=$(echo $line | awk -v OFS="\t" '{print $4}')
    ID=~{name}
    CHR2=$(echo $line | awk -v OFS="\t" '{print $3}')
    # echo "Tloc:" $ID $BP1_interval "and" $BP2_interval
    zcat ~{tloc_pe_evidence} | awk -v OFS="\t" -v var2=$CHR2 -v pos1=$(($BP2-1000)) -v pos2=$(($BP2+1000)) '{if ($4 == var2 && $5>pos1 && $5<pos2) print$0}' | awk -v OFS="\t" '{print $5}' | sort -n -k 1,1 > PE_list
    zcat ~{tloc_pe_evidence} | awk -v OFS="\t" -v var2=$CHR2 -v pos1=$(($BP2-1000)) -v pos2=$(($BP2+1000)) '{if ($4 == var2 && $5>pos1 && $5<pos2) print$0}' > $ID.list
    start=$(awk 'NR==1' PE_list)
    end=$(awk 'END{print}' PE_list)
    cat $ID.list | awk -v OFS="\t" '{print $7}' | sort | uniq -c | sort -nr > $ID.size
    Sample_size=$(cat $ID.list | awk -v OFS="\t" '{print $7}' | sort | uniq -c | sort -nr | wc -l)
    Sample_size_min4=$(cat $ID.list | awk -v OFS="\t" '{print $7}' | sort | uniq -c | sort -nr | awk '{if ($1 > 3) print $0}' | wc -l)
    Sample_size_min10=$(cat $ID.list | awk -v OFS="\t" '{print $7}' | sort | uniq -c | sort -nr | awk '{if ($1 > 9) print $0}' | wc -l)
    # echo "Translocation $ID distance from $CHR2 : $BP2 starts at $start ends at $end" in $Sample_size > $ID.report
    cat $ID.list $ID.size > ~{name}.info
    rm $ID.size
    rm $ID.list
    # rm $ID.report
    awk -v line=$i -v OFS="\t" '{print $0}' ~{tloc_bed} | sed "s/$/\t$Sample_size/" | sed "s/$/\t$Sample_size_min4/" | sed "s/$/\t$Sample_size_min10/" > ~{name}.freq
    done

    rm PE_list

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

task merge_frequencies {
  input {
    Array[File] tloc_freq_files
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
    File tloc_merged_freq ="cohort.freq.gz"
  }

  command <<<
    set -ex
    cat ~{sep=" " tloc_freq_files} | bgzip > cohort.freq.gz
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

task create_summary {
  input {
    File tloc_sample_bed
    Array[File] tloc_info_file
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
    File tloc_summary ="~{name}.summary"
  }

  command <<<
    set -ex

    grep ~{name} ~{tloc_sample_bed} > ~{name}_samples.bed

    caselist=$1
    i=0
    cat ~{tloc_sample_bed} | while read line
    do
      i=$((i+1))
      sample_ID=$(echo $line | awk -v OFS="\t" '{print $6}')
      tloc_ID=$(echo $line | awk -v OFS="\t" '{print $5}')
      BP1_1_plus=$(cat ~{tloc_info_file} | grep -v $sample_ID | awk '{if($3 == "+") print $2}' | sort  | awk 'NR==1 {print}')
      BP1_2_plus=$(cat ~{tloc_info_file} | grep -v $sample_ID | awk '{if($3 == "+") print $2}' | sort  | awk 'END {print}')
      BP1_1_minus=$(cat ~{tloc_info_file} | grep -v $sample_ID | awk '{if($3 == "-") print $2}' | sort  | awk 'NR==1 {print}')
      BP1_2_minus=$(cat ~{tloc_info_file} | grep -v $sample_ID | awk '{if($3 == "-") print $2}' | sort  | awk 'END {print}')
      BP1_int_plus=$((BP1_1_plus-BP1_2_plus))
      BP1_int_minus=$((BP1_1_minus-BP1_2_minus))
      BP2_1_plus=$(cat ~{tloc_info_file} | grep -v $sample_ID | awk '{if($6 == "+") print $5}' | sort  | awk 'NR==1 {print}')
      BP2_2_plus=$(cat ~{tloc_info_file} | grep -v $sample_ID | awk '{if($6 == "+") print $5}' | sort  | awk 'END {print}')
      BP2_1_minus=$(cat ~{tloc_info_file} | grep -v $sample_ID | awk '{if($6 == "-") print $5}' | sort  | awk 'NR==1 {print}')
      BP2_2_minus=$(cat ~{tloc_info_file} | grep -v $sample_ID | awk '{if($6 == "-") print $5}' | sort  | awk 'END {print}')
      BP2_int_plus=$((BP2_1_plus-BP2_2_plus))
      BP2_int_minus=$((BP2_1_minus-BP2_2_minus))
      awk -v line=$i -v OFS="\t" '{print $0}' ~{tloc_sample_bed} | \
        sed "s/$/\t$BP1_int_plus/" | sed "s/$/\t$BP1_int_minus/" | \
        sed "s/$/\t$BP2_int_plus/" | sed "s/$/\t$BP2_int_minus/" | \
        sed "s/$/\t$BP1_1_plus/" | sed "s/$/\t$BP1_2_plus/" | \
        sed "s/$/\t$BP1_1_minus/" | sed "s/$/\t$BP1_2_minus/" |  \
        sed "s/$/\t$BP2_1_plus/" | sed "s/$/\t$BP2_2_plus/" | \
        sed "s/$/\t$BP2_1_minus/" | sed "s/$/\t$BP2_2_minus/" > ~{name}.summary
    done

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

task merge_summary {
  input {
    Array[File] tloc_summary_files
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
    File tloc_merged_summary ="cohort.summary.gz"
  }

  command <<<
    set -ex
    cat ~{sep=" " tloc_summary_files} | bgzip > cohort.summary.gz
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