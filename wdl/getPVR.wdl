version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow getPVR {
  input {
    File vcf_metrics_tsv
    File sample_crai_cram
    File trio_uri
    String pvr_docker
  }

  call getFamList {
    input:
      trio_uri=trio_uri,
      pvr_docker=pvr_docker
  }

  scatter (fam_id in getFamList.fam_ids) {
    call subsetFam {
      input:
        vcf_metrics_tsv=vcf_metrics_tsv,
        trio_uri=trio_uri,
        sample_crai_cram=sample_crai_cram,
        fam_id=fam_id,
        pvr_docker=pvr_docker
    }
    call calculatePVR {
      input:
        indel_file=subsetFam.indel_file,
        trio_uri=trio_uri,
        fam_crais=subsetFam.fam_crais,
        fam_crams=subsetFam.fam_crams,
        pvr_docker=pvr_docker
    }
  }  

  call combinePVR {
    input:
      vcf_metrics_tsv=vcf_metrics_tsv,
      getPVR_out=calculatePVR.getPVR_out,
      getPVR_redo=calculatePVR.redo,
      pvr_docker=pvr_docker
  }  

  output {
    Array[File] getPVR_out = calculatePVR.getPVR_out
    Array[File] getPVR_redo = calculatePVR.redo
    File final_getPVR_out = combinePVR.final_getPVR_out
    File final_getPVR_redo = combinePVR.final_getPVR_redo
  }
}

task getFamList {
  input {
    File trio_uri
    String pvr_docker
  }

  runtime {
    docker: pvr_docker
  }

  command {
    cat ~{trio_uri} | tail -n +2 | cut -f1 | uniq > fam_list.txt
  }

  output {
    Array[String] fam_ids = read_lines("fam_list.txt")
  }
}

task subsetFam {
  input {
    File vcf_metrics_tsv
    File trio_uri
    File sample_crai_cram
    String fam_id
    String pvr_docker
  }

  runtime {
    docker: pvr_docker
  }

  command {
    head -1 ~{vcf_metrics_tsv} > ~{fam_id}.txt
    awk -v fam_id=~{fam_id} '$10==fam_id' ~{vcf_metrics_tsv} >> ~{fam_id}.txt
    awk -v fam_id=~{fam_id} '$1==fam_id' ~{trio_uri} | cut -f3-5 | tr '\t' '\n' | uniq > samples.txt
    grep -f samples.txt ~{sample_crai_cram} | cut -f2 > crai_files.txt
    grep -f samples.txt ~{sample_crai_cram} | cut -f3 > cram_files.txt
  }

  output {
    File indel_file = fam_id + '.txt'
    Array[File] fam_crais = read_lines("crai_files.txt")
    Array[File] fam_crams = read_lines("cram_files.txt")
  }
}

task calculatePVR {
  input {
    File indel_file
    File trio_uri
    Array[File] fam_crams
    Array[File] fam_crais
    String pvr_docker
    RuntimeAttr? runtime_attr_override
  }

  # Runtime parameters adapted from gatk-sv "CollectCoverage.wdl"
  Int num_cpu = 4
  Int mem_size_gb = 16
  Int vm_disk_size = 300

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    docker: pvr_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  String fam_id = basename(indel_file, ".txt")
  String pvr_filename = "${fam_id}_with_PVR.txt"
  String pvr_redo = "${fam_id}_redo.txt"

  command {
    set -euo pipefail
    
    perl /home/get_pvr_from_crams_wdl.pl -i ~{indel_file} -m ~{trio_uri} -b ${sep="\\;" fam_crams}

  }
  output {
    File getPVR_out = pvr_filename
    File redo = pvr_redo
  }
}

task combinePVR {
  input {
    File vcf_metrics_tsv
    Array[File] getPVR_out
    Array[File] getPVR_redo 
    String pvr_docker
    RuntimeAttr? runtime_attr_override
  }

  # Runtime parameters adapted from gatk-sv "CollectCoverage.wdl"
  Int num_cpu = 4
  Int mem_size_gb = 16
  Int vm_disk_size = 300

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    docker: pvr_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  String getPVR_out_filename = basename(vcf_metrics_tsv, '.tsv') + '_with_PVR.tsv'
  String getPVR_redo_filename = basename(vcf_metrics_tsv, '.tsv') + '_PVR_redo.tsv'
  command {
    cat ~{vcf_metrics_tsv} | head -1 > ~{getPVR_redo_filename}
    echo -e $(cat ~{getPVR_redo_filename})"\tChild_tag\tall_count\tref_count\talt_count ref_prop\talt_prop\tsam_real_alt\treal_alt_count\treal_alt_prop complex_count Father_tag\tfa_all_count\tfa_ref_count\tfa_alt_count\tfa_ref_prop fa_alt_prop fa_sam_real_alt fa_real_alt_count fa_real_alt_prop\tfa_complex_count\tMother_tag\tmo_all_count\tmo_ref_count\tmo_alt_count\tmo_ref_prop mo_alt_prop mo_sam_real_alt mo_real_alt_count mo_real_alt_prop\tmo_complex_count" > ~{getPVR_out_filename} 

    cat ~{sep=' ' getPVR_out} >> ~{getPVR_out_filename}    
    cat ~{sep=' ' getPVR_redo} >> ~{getPVR_redo_filename}
  }

  output {
    File final_getPVR_out = getPVR_out_filename
    File final_getPVR_redo = getPVR_redo_filename
  }
}