version 1.0

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
workflow pathway_enrichment {
  input {
    # File genes_file
    File info_file #Fu et al. 2022 gene list with HGNC symbols and Entrez IDs
    File pedigree_file
    File denovo_file
    File phenotype_pairs
    File hgnc_file
    Int eigenvalue
    String mutation_type
    String docker_pathway_enrichment
    RuntimeAttr? runtime_attr_pathway_enrichment_00
    RuntimeAttr? runtime_attr_pathway_enrichment_01
    RuntimeAttr? runtime_attr_pathway_enrichment_02
    RuntimeAttr? runtime_attr_pathway_enrichment_merge
  }

  ## Read the tab-separated file as an Array of Arrays
  Array[Array[String]] pairs = read_tsv(phenotype_pairs)

  call pathway_enrichment_00 {
      input:
        info_file = info_file,
        hgnc_file = hgnc_file,
        docker_path = docker_pathway_enrichment,
        runtime_attr_override = runtime_attr_pathway_enrichment_00
    }

  scatter (row in pairs) {
    # Access the first and second columns directly via array indexing
    String phenotype1_raw = row[0]
    String phenotype2_raw = row[1]

    # Sanitize by replacing / and space with _
    String phenotype1 = sub(sub(phenotype1_raw, "/", "_"), " ", "_")
    String phenotype2 = sub(sub(phenotype2_raw, "/", "_"), " ", "_")

    # Create the outfile name
    # String outfile = "${phenotype1}_${phenotype2}.RData"

    call pathway_enrichment_01 {
      input:
        info_file = pathway_enrichment_00.genelist_entrez_id,
        pedigree_file = pedigree_file,
        denovo_file = denovo_file,
        phenotype1 = phenotype1_raw,
        phenotype2 = phenotype2_raw,
        docker_path = docker_pathway_enrichment,
        runtime_attr_override = runtime_attr_pathway_enrichment_01
    }

    call pathway_enrichment_02 {
      input:
        marker_file = pathway_enrichment_01.marker_counts, # Output from pathway_enrichment_01
        genes_file = pathway_enrichment_01.dnv_counts,
        phenotype1_raw = phenotype1_raw,
        phenotype2_raw = phenotype2_raw,
        phenotype1 = phenotype1,
        phenotype2 = phenotype2,
        mutation_type = mutation_type,
        eigenvalue = eigenvalue,
        docker_path = docker_pathway_enrichment,
        runtime_attr_override = runtime_attr_pathway_enrichment_02        
    }
  }

  call merge_pathway_enrichment_outputs {
      input:
        go_pathway_output = pathway_enrichment_02.go_pathway_output,
        go_volcano_plot = pathway_enrichment_02.go_volcano_plot,
        docker_path = docker_pathway_enrichment,
        runtime_attr_override = runtime_attr_pathway_enrichment_merge        

    }

  output {
    File merged_outputs =  merge_pathway_enrichment_outputs.merged_output
  }
}

# TASK DEFINITIONS

task pathway_enrichment_00 {
  input {
    File info_file
    File hgnc_file
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
    File genelist_entrez_id = "genelist_fixed_entrez_id.txt"
  }

  command <<<
    set -ex

    Rscript /src/variant-interpretation/scripts/pathway_enrichment-00.R \
      --genes "~{info_file}" \
      --id_hgnc "~{hgnc_file}"

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

task pathway_enrichment_01 {
  input {
    File info_file
    File pedigree_file
    File denovo_file
    String phenotype1
    String phenotype2
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
    File marker_counts = "hpo_marker_pair_counts.txt"
    File dnv_counts = "crossdev_pair_dn_counts.txt"
  }

  command <<<
    set -ex

    Rscript /src/variant-interpretation/scripts/pathway_enrichment-01.R \
      --info "~{info_file}" \
      --pedigree "~{pedigree_file}" \
      --denovo "~{denovo_file}"  \
      --phenotype1 "~{phenotype1}" \
      --phenotype2 "~{phenotype2}"
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

task pathway_enrichment_02 {
  input {
    File marker_file
    File genes_file
    Int eigenvalue
    String phenotype1_raw
    String phenotype2_raw
    String phenotype1
    String phenotype2
    String mutation_type
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
    File go_pathway_output = "GO_pathway_~{phenotype1}_~{phenotype2}_~{mutation_type}.tsv"
    File go_volcano_plot = "Volcano_~{phenotype1}_~{phenotype2}_~{mutation_type}.pdf"
  }

  command <<<
    set -ex

    Rscript /src/variant-interpretation/scripts/pathway_enrichment-02.R \
      --phenotype1 "~{phenotype1_raw}" \
      --phenotype2 "~{phenotype2_raw}" \
      --mutation "~{mutation_type}" \
      --marker "~{marker_file}" \
      --genes "~{genes_file}" \
      --eigenvalue "~{eigenvalue}"
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

task merge_pathway_enrichment_outputs {
  input {
    Array[File] go_pathway_output
    Array[File] go_volcano_plot
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
    File merged_output = "pathway_enrichment_results.tar.gz"
  }

  command <<<
    set -ex

    # Create a directory for the outputs if it doesn't exist
    mkdir -p pathway_enrichment_results

    # Move the output files to the directory with a standardized naming convention
    mv ~{sep=" " go_pathway_output} pathway_enrichment_results/
    mv ~{sep=" " go_volcano_plot} pathway_enrichment_results/

    tar -czf pathway_enrichment_results.tar.gz pathway_enrichment_results
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