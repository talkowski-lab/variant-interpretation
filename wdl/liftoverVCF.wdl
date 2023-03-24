version 1.0

## This workflow performs lifover of a VCF to a new genome assembly
## Author: asanchis@broadinstitute.org

# WORKFLOW DEFINITION
workflow liftoverVCF {
    input {
        File input_vcf
        File chain_file
        File new_ref_genome
        String output_name
        String docker_path
    }

    call liftover {
        input:
            input_vcf = input_vcf,
            chain_file = chain_file,
            new_ref_genome = new_ref_genome,
            output_name = output_name,
            docker = docker_path
    }

    output {
        File output_vcf = liftover.output_name
        File output_vcf_index = liftover.output_vcf_index
    }
}

# TASK DEFINITIONS

task liftover {
    input {
        File input_vcf
        File output_name
        String chain_file
        String new_ref_genome
        String docker_path

        # Runtime parameters
        Int? mem_gb
        Int? disk_space_gb
        Int? preemptible_attempts
    }

    Boolean use_ssd = false
    Int machine_mem_gb = select_first([mem_gb, 3])
    Int command_mem_gb = machine_mem_gb - 1

    command {
        set -e

        gatk LiftoverVcf \
            -I ~{input_vcf} \
            -O ~{output_name} \
            --CHAIN ~{chain_file} \
            --REJECT rejected.vcf.gz \
            -R ~{genome_build}

        tabix -p vcf ~{output_name}

  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + select_first([disk_space_gb, 20]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
  }
  output {
    File output_name = output_name
    File output_name_index = "{output_name}.tbi"
  }
}