version 1.0

## This workflow reformat the output of GATK-SV to match seqr requirements

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
workflow ReformatGregor {
  input {
    File concordance_vcf
    File annotated_vcf
    File? annotated_vcf_idx
    File pedigree
    File non_PAR
    File new_ids
    String prefix
    String docker_gregor
    RuntimeAttr? runtime_attr_gregor
    RuntimeAttr? runtime_attr_concordance_ids
  }

  call fix_chrX_GT{
      input:
        vcf = concordance_vcf,
        non_PAR = non_PAR,
        prefix = prefix,
        pedigree = pedigree,
        docker_path = docker_gregor,
        runtime_attr_override = runtime_attr_gregor
  }

  call change_sampleIDs{
      input:
        vcf = fix_chrX_GT.out_vcf,
        new_ids = new_ids,
        docker_path = docker_gregor,
        prefix = prefix,
        runtime_attr_override = runtime_attr_gregor
  }

  call fixGnomAD{
      input:
        vcf = change_sampleIDs.out_vcf,
        prefix = prefix,
        docker_path = docker_gregor,
        runtime_attr_override = runtime_attr_gregor
  }

  call strvctvre{
      input:
        vcf = fixGnomAD.out_vcf,
        prefix = prefix,
        docker_path = docker_gregor,
        runtime_attr_override = runtime_attr_gregor
  }

  call fixConcordanceIDs{
      input:
        vcf = strvctvre.out_vcf,
        prefix = prefix,
        docker_path = docker_gregor,
        runtime_attr_override = runtime_attr_concordance_ids
  }

  call annotateFilter{
      input:
        annotated_vcf = annotated_vcf,
        annotated_vcf_idx = annotated_vcf_idx,
        vcf = fixConcordanceIDs.out_vcf,
        vcf_index = fixConcordanceIDs.out_index,
        prefix = prefix,
        docker_path = docker_gregor,
        runtime_attr_override = runtime_attr_gregor
  }

  call cleanHeader{
      input:
        vcf = annotateFilter.out_vcf,
        prefix = prefix,
        docker_path = docker_gregor,
        runtime_attr_override = runtime_attr_gregor
  }

  output {
    File seqr_vcf = cleanHeader.out_vcf
    File seqr_vcf_index = cleanHeader.out_index
    }

}

# TASK DEFINITIONS
task fix_chrX_GT {
  input {
    File vcf
    File non_PAR
    File pedigree
    File docker_path
    String prefix
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
    File out_vcf = "~{prefix}.fixChrX.vcf.gz"
  }

  command {
    set -e
    bedtools intersect -a ~{vcf} -b ~{non_PAR} -f 0.5 | cut -f3 > chrX_ids.txt
    python3 fix_GT_chrXmale.py --vcf ~{vcf} --out ~{prefix}.fixChrX.vcf.gz --ped ~{pedigree} --ids chrX_ids.txt
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

task change_sampleIDs {
  input {
    File vcf
    File new_ids
    File docker_path
    String prefix
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
    File out_vcf = "~{prefix}.sampleIDs.vcf.gz"
  }

  command {
    set -e
    bcftools reheader ~{vcf} -s ~{new_ids} > ~{prefix}.sampleIDs.vcf.gz
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

task fixGnomAD {
  input {
    File vcf
    File docker_path
    String prefix
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
    File out_vcf = "~{prefix}.fixHeader.vcf.gz"
  }

  command {
    set -e
    bcftools view ~{vcf} | sed -e 's/gnomAD_V2_AC_AF/gnomAD_V2_AC/g' | \
      sed -e 's/gnomAD_V2_AN_AF/gnomAD_V2_AN/g' | \
      bcftools view -O z -o ~{prefix}.fixHeader.vcf.gz
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

task strvctvre {
  input {
    File vcf
    File docker_path
    String prefix
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
    File out_vcf = "~{prefix}.StrVCTVRE.vcf.gz"
  }

  command {
    set -e
    cd ~/StrVCTVRE
    python3 StrVCTVRE.py -i ~{vcf} -o ~{prefix}.StrVCTVRE.vcf.gz
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

task fixConcordanceIDs {
  input {
    File vcf
    File docker_path
    String prefix
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
    File out_vcf = "~{prefix}.concorID.vcf.gz"
    File out_index = "~{prefix}.concorID.vcf.gz.tbi"
  }

  command {
    set -e
    python3 updateConcordanceID.py --input ~{vcf} --output ~{prefix}.concorID.vcf
    bgzip ~{prefix}.concorID.vcf
    tabix -p vcf ~{prefix}.concorID.vcf
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

task annotateFilter {
  input {
    File annotated_vcf
    File? annotated_vcf_idx
    File vcf
    File vcf_index
    File docker_path
    String prefix
    RuntimeAttr? runtime_attr_override
  }

  File annotated_vcf_idx = select_first([annotated_vcf_idx, annotated_vcf + ".tbi"])


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
    File out_vcf = "~{prefix}.concorID.vcf.gz"
    File out_index = "~{prefix}.concorID.vcf.gz.tbi"
  }

  command {
    set -e
    bcftools annotate -c FILTER -a ~{annotated_vcf} -O z -o ~{prefix}.fixFilter.vcf.gz ~{vcf}
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

task cleanHeader {
  input {
    File vcf
    File docker_path
    String prefix
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
    File out_vcf = "~{prefix}.seqr.vcf.gz"
    File out_index = "~{prefix}.seqr.vcf.gz.tbi"
  }

  command {
    set -e

    ##Get VCF header
    bcftools view -h ~{vcf} > header.txt

    ##Remove bcftools view commands
    grep -v bcftools_viewCommand header.txt > header_fix.txt

    ##Reheader
    bcftools reheader -h header_fix.txt -o ~{prefix}.seqr.vcf.gz ~{vcf}

    ##Make index
    tabix -p vcf ~{prefix}.seqr.vcf.gz
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