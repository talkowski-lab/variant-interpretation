## This script extract GD overlap from VCF and Raw depth data - designed for GREGoR cohort

version 1.0

# IMPORT STRUCTS
struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow GenomicDisorders {
  input {
    File vcf
    File genomic_disorders
    File raw_depth
    File raw_files_list
    Array[String] contigs
    String prefix

    String docker_genomic_disorders

    RuntimeAttr? runtime_attr_gd
  }

  Array[String] raw_files = transpose(read_tsv(raw_files_list))[1]

  call reformatVCF{
      input:
        vcf = vcf,
        prefix = prefix,
        docker_path = docker_genomic_disorders,
        runtime_attr_override = runtime_attr_gd
  }

  call getVCFoverlap{
    input:
      bed = reformatVCF.out_bed,
      genomic_disorders = genomic_disorders,
      prefix = prefix,
      docker_path = docker_genomic_disorders,
      runtime_attr_override = runtime_attr_gd
  }

  scatter(raw_file in raw_files){
    call raw_reformatVCF {
      input:
        vcf = raw_file,
        docker_path = docker_genomic_disorders,
        prefix = basename(raw_file, ".vcf.gz"),
        runtime_attr_override = runtime_attr_gd
        }
    }

    call raw_mergeBed {
      input:
      bed_files = raw_reformatVCF.bed_output,
      docker_path = docker_genomic_disorders,
      runtime_attr_override = runtime_attr_gd
    }

    scatter (contig in contigs){
      call raw_divideByChrom {
        input:
          bed_file = raw_mergeBed.concat_bed_output,
          chromosome = contig,
          docker_path = docker_genomic_disorders,
          runtime_attr_override = runtime_attr_gd
      }

      call raw_reformatBedDepth{
        input:
          per_chromosome_bed_file = raw_divideByChrom.per_chromosome_bed_output,
          chromosome=contig,
          docker_path = docker_genomic_disorders,
          runtime_attr_override = runtime_attr_gd
      }
    }

  output {
    File vcf_to_bed = reformatVCF.out_bed
    File raw_reformatted = raw_reformatBedDepth.reformatted_output
    }

}

# TASK DEFINITIONS
task reformatVCF {
  input {
    File vcf
    String docker_path
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
    File out_bed = "~{prefix}.bed.gz"
    File out_ref_bed = "{prefix}.ref.bed.gz"
  }

  command <<<
    set -e
    svtk vcf2bed -i ALL --include-filters ~{vcf} - | bgzip -c > ~{prefix}.bed.gz
    zcat ~{prefix}.bed.gz | \
      grep -E "DEL|DUP" | \
      awk '{print $1"_"$5"\t"$2"\t"$3"\t"$4"\t"$5}' | \
      grep -v ^#chrom | \
      bgzip -c > ~{prefix}.ref.bed.gz

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

task getVCFoverlap {
  input {
    File bed
    File genomic_disorders
    String docker_path
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
    File out_bed = "~{bed}.gd.bed.gz"
  }

  command <<<
    set -e
    bedtools intersect -wa -wb -f 0.3 -r -a ~{bed} -b ~{genomic_disorders} | bgzip -c > ~{bed}.gd.bed.gz

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

task raw_reformatVCF{
    input{
        File vcf
        String docker_path
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + input_size * 2.0),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File bed_output = "${prefix}.bed.gz"
    }

    command {
        set -euo pipefail

        #convert from vcf to bed file
        svtk vcf2bed ~{vcf_file} --info SVTYPE - | bgzip -c > ${prefix}.bed.gz
    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker_path
    }
}

task raw_mergeBed{
    input{
        Array[File] bed_files
        String docker_path
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(bed_files, "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + input_size * 1.5),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File concat_bed_output = "concat.bed.gz"
    }

    command {
        set -euo pipefail
        zcat ${sep=" " bed_files} | bgzip -c > concat.bed.gz

    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker_path
    }
}

task raw_divideByChrom{
    input{
        File bed_file
        String chromosome
        String docker_path
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(bed_file, "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + input_size * 1.3),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File per_chromosome_bed_output = "${chromosome}.bed.gz"
    }

    command {
        set -euo pipefail

        zcat ${bed_file} | \
        grep -w ^${chromosome} | \
        bgzip -c > ${chromosome}.bed.gz

    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker_path
    }
}

task raw_reformatBedDepth{
    input{
        File per_chromosome_bed_file
        String chromosome
        String docker_path
        RuntimeAttr? runtime_attr_override
    }

    Float bed_file_size = size(per_chromosome_bed_file, "GB")
    Float ped_size = size(ped_input, "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + ped_size + bed_file_size * 2.0),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
       File reformatted_output = "${chromosome}.reformatted.sorted.bed.gz"
    }

    command {
        set -euo pipefail

        #reformat bed file
        Rscript /src/variant-interpretation/scripts/reformatRawBed_all.R ${per_chromosome_bed_file} ${chromosome}.reformatted.bed
        sortBed -i ${chromosome}.reformatted.bed | bgzip -c > ${chromosome}.reformatted.sorted.bed.gz
    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker_path
    }
}