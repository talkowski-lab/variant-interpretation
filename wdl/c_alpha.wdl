## Running C-Alpha

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
workflow calpha {
  input {
    File phenotype_pairs
    File pedigree_file
    File snvs_indels
    File genes_file
    String docker_calpha
    RuntimeAttr? runtime_attr_calpha
  }

    ##Here, read the file phenotype_pairs, then scatter
    Array[String] pairs = read_lines(phenotype_pairs)

    scatter (pair in pairs) {
      # Split the tab-delimited pair into two phenotypes
      Array[String] split_pair = split(pair, "\t")
      String phenotype1_raw = split_pair[0]
      String phenotype2_raw = split_pair[1]

      # Sanitize by replacing / and space with _
      String phenotype1 = sub(sub(phenotype1_raw, "/", "_"), " ", "_")
      String phenotype2 = sub(sub(phenotype2_raw, "/", "_"), " ", "_")

      # Create the outfile name
      String outfile = "${phenotype1}_${phenotype2}.RData"

      call calpha_task {
        input:
          phenotype1 = phenotype1_raw,
          phenotype2 = phenotype2_raw,
          pedigree_file = pedigree_file,
          snvs_indels = snvs_indels,
          genes_file = genes_file,
          outfile = outfile,
          docker_path = docker_calpha,
          runtime_attr_override = runtime_attr_calpha
  }
}

  output {
    Array[File] output_calpha = calpha_task.rdata_output
  }
}

# TASK DEFINITIONS

task calpha_task {
  input {
    String phenotype1
    String phenotype2
    String outfile
    File pedigree_file
    File snvs_indels
    File genes_file
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
    File rdata_output = "~{outfile}"
  }

  command <<<
    set -ex

    Rscript /src/variant-interpretation/scripts/c_alpha.R \
      --phenotype1 ~{phenotype1} \
      --phenotype2 ~{phenotype2} \
      --pedigree ~{pedigree_file} \
      --input_snv ~{snvs_indels} \
      --genes ~{genes_file} \
      --output ~{outfile}
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