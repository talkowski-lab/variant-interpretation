version 1.0

workflow ReblockGVCFs {

  String pipeline_version = "2.0.3"

  input {
    File gvcf
    File gvcf_index
    File collaborator_participant_id
    File? calling_interval_list
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Float? tree_score_cutoff
    String? annotations_to_keep_command
    String? annotations_to_remove_command
    Boolean? move_filters_to_genotypes
    Boolean move_bam_or_cram_files=false
    Boolean disable_sequence_dictionary_validation=false
    String docker_image
    Array[String] exclude_contigs=[]
    String gvcf_file_extension = ".g.vcf.gz"
  }

  String gvcf_basename = basename(gvcf, gvcf_file_extension)
  
  
  call LocalizeReads {
    input:
      reads_path = gvcf,
      reads_index = gvcf_index,
      collaborator_participant_id = collaborator_participant_id,
      move_files = move_bam_or_cram_files,
      docker_image = docker_image
    }

  if (length(exclude_contigs)!=0) {
    call removeExtraContigs {
        input:
        gvcf=LocalizeReads.output_file,
        gvcf_index=LocalizeReads.output_index,
        docker_image=docker_image,
        exclude_contigs=exclude_contigs
    }
  }

  call Reblock {
    input:
      gvcf = select_first([removeExtraContigs.output_vcf, LocalizeReads.output_file]),
      gvcf_index = select_first([removeExtraContigs.output_vcf_index, LocalizeReads.output_index]),
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      tree_score_cutoff = tree_score_cutoff,
      annotations_to_keep_command = annotations_to_keep_command,
      disable_sequence_dictionary_validation=disable_sequence_dictionary_validation,
      output_vcf_filename = gvcf_basename + ".rb.g.vcf.gz"
  }


  output {
    File output_vcf = Reblock.output_vcf
    File output_vcf_index = Reblock.output_vcf_index
  }
  meta {
    allowNestedInputs: true
  }
}

task LocalizeReads {
  input {
    File reads_path
    File reads_index
    String collaborator_participant_id
    Boolean move_files = false
    String docker_image = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
  }

  Float input_size = if move_files then size(reads_path, "GB") * 2 else size(reads_path, "GB")
  
  runtime {
    memory: "3750 MiB"
    disks: "local-disk " + disk_size + " HDD"
    bootDiskSizeGb: 15
    preemptible: 3
    docker: docker_image
  }

  Int disk_size = ceil(50 + size(reads_path, "GB"))

  command {
    set -exuo pipefail

    # When this pipeline is run on an HPC, moving files could lead to
    # moving the files from their original source, compared to moving
    # them from one directory of the VM to another when run on Cloud.
    # Therefore, to avoid moving files unexpectedly, we provide both
    # options for moving and copying, and set the copy as default.
    # Note that, when copying the files, the task can be slower depending
    # on the file size and IO performance and will need additional disk
    # space, hence it will be more expensive to run.

    if ~{move_files}; then
      mv ~{reads_path} $(basename ~{reads_path})
      mv ~{reads_index} $(basename ~{reads_index})
    else
      cp ~{reads_path} $(basename ~{reads_path})
      cp ~{reads_index} $(basename ~{reads_index})
    fi
  }
  output {
    File output_file = basename(reads_path)
    File output_index = basename(reads_index)
  }
}

task removeExtraContigs {
  input {
    File gvcf
    File gvcf_index
    String docker_image = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    Int additional_disk = 20
    Array[String] exclude_contigs
  }

  Int disk_size = ceil((size(gvcf, "GiB")) * 4) + additional_disk
  String file_ext = if sub(basename(gvcf), '\\.g\\.vcf\\.gz', '')!=basename(gvcf) then '.gz.vcf.gz' else '.gvcf.gz'
  String output_filename = basename(gvcf, file_ext) + '.rm.extra.contigs.vcf.gz'

  command {
    set -eou pipefail
    bcftools index -s ~{gvcf} | grep -Fv ~{sep=' -e ' exclude_contigs} > contigs.txt
    bcftools view -r $(cat contigs.txt | cut -f1 | tr '\n' ',') -Oz -o ~{output_filename} ~{gvcf}
    bcftools index -t ~{output_filename}
  }

  runtime {
    memory: "3750 MiB"
    disks: "local-disk " + disk_size + " HDD"
    bootDiskSizeGb: 15
    preemptible: 3
    docker: docker_image
  }

  output {
    File output_vcf = output_filename
    File output_vcf_index = output_filename + '.tbi'
  }
}

task Reblock {

  input {
    File gvcf
    File gvcf_index
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    String output_vcf_filename
    String docker_image = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    Int additional_disk = 20
    Boolean disable_sequence_dictionary_validation
    String? annotations_to_keep_command
    Float? tree_score_cutoff
  }

  Int disk_size = ceil((size(gvcf, "GiB")) * 4) + additional_disk

  command {
    set -eou pipefail

    gatk --java-options "-Xms3000m -Xmx3000m" \
      ReblockGVCF \
      -R ~{ref_fasta} \
      -V ~{gvcf} \
      -do-qual-approx \
      --floor-blocks -GQB 20 -GQB 30 -GQB 40 \
      ~{annotations_to_keep_command} \
      ~{if disable_sequence_dictionary_validation then "--disable-sequence-dictionary-validation true" else ""} \
      ~{"--tree-score-threshold-to-no-call " + tree_score_cutoff} \
      -O ~{output_vcf_filename}
  }

  runtime {
    memory: "3750 MiB"
    disks: "local-disk " + disk_size + " HDD"
    bootDiskSizeGb: 15
    preemptible: 3
    docker: docker_image
  }

  output {
    File output_vcf = output_vcf_filename
    File output_vcf_index = output_vcf_filename + ".tbi"
  }
}
