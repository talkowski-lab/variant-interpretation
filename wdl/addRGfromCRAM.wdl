version 1.0

## Copyright Broad Institute, 2019
##
## Adaptation of multiple worflows to add RG to input CRAM
## Author: asanchis@broadinstitute.org
##
## Adapted from:
## CRAM to BAM: github.com/gatk-workflows/seq-format-conversion/CRAM-to-BAM:master
## addOrReplaceReadGroups: margolis/addOrReplaceReadGroups/1
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.

# WORKFLOW DEFINITION

workflow AddOrReplaceReadGroupsWorkflow {
    input {
        File ref_fasta
        File ref_fasta_index_fai
        File ref_fasta_index_gzi
        File ref_dict
        File bam_or_cram
        String sample_name
        String gotc_docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
        Int preemptible_tries = 3
        String readgroup_library
        String readgroup_platform
        String readgroup_run_barcode
        Int memory_addOrReplaceGroups
        Int disk_addOrReplaceGroups
    }

    Boolean is_cram =
        basename(bam_or_cram, ".cram") + ".cram" == basename(bam_or_cram)

    if (is_cram){
        #Converts CRAM to SAM to BAM and makes BAI
        call CramToBamTask {
            input:
                ref_fasta = ref_fasta,
                ref_fasta_index_fai = ref_fasta_index_fai,
                ref_fasta_index_gzi = ref_fasta_index_gzi,
                ref_dict = ref_dict,
                bam_or_cram = bam_or_cram,
                sample_name = sample_name,
                docker_image = gotc_docker,
                preemptible_tries = preemptible_tries
        }
    }

    File bam_file = select_first([CramToBamTask.outputBam, bam_or_cram])

    ##Adds or replaces read groups and makes BAI
    call AddOrReplaceReadGroups {
        input:
            inputBam = bam_file,
            sample_name = sample_name,
            readgroupLibrary = readgroup_library,
            readgroupPlatform = readgroup_platform,
            readgroupRunBarcode = readgroup_run_barcode,
            memoryGb = memory_addOrReplaceGroups,
            diskSpaceGb = disk_addOrReplaceGroups
    }
    #Validates Bam after add or replace groups
    call ValidateSamFile {
        input:
            input_bam = AddOrReplaceReadGroups.bamWithReadGroupAdded,
            docker_image = gotc_docker,
            preemptible_tries = preemptible_tries
    }
    # Outputs that will be retained when execution is complete
    output {
    	File bam_with_RG = AddOrReplaceReadGroups.bamWithReadGroupAdded
        File bai_with_RG = AddOrReplaceReadGroups.bamWithReadGroupAddedIndex
        File validate_sam_report = ValidateSamFile.report
    }
}

#Task Definitions

task CramToBamTask {
    input {
        # Command parameters
        File ref_fasta
        File ref_fasta_index_fai
        File ref_fasta_index_gzi
        File ref_dict
        File bam_or_cram
        String sample_name

        # Runtime parameters
        Int addtional_disk_size = 20
        Int machine_mem_size = 15
        String docker_image
        Int preemptible_tries
    }
    Float output_bam_size = size(bam_or_cram, "GB") / 0.60
    Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index_fai, "GB") + size(ref_dict, "GB")
    Int disk_size = ceil(size(bam_or_cram, "GB") + output_bam_size + ref_size) + addtional_disk_size

    #Calls samtools view to do the conversion
    command {
        set -eo pipefail

        samtools view -h -T ~{ref_fasta} ~{bam_or_cram} | \
        samtools view -b -o ~{sample_name}.bam -
    }

    #Run time attributes:
    #Use a docker with samtools. Set this up as a workspace attribute.
    #cpu of one because no multi-threading is required. This is also default, so don't need to specify.
    #disk_size should equal input size + output size + buffer
    runtime {
        docker: docker_image
        memory: machine_mem_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible_tries
    }

    #Outputs a BAM and BAI with the same sample name
    output {
        File outputBam = "~{sample_name}.bam"
    }
}


#Validates BAM output to ensure it wasn't corrupted during the file conversion

task ValidateSamFile {
    input {
        File input_bam
        Int addtional_disk_size = 10
        Int machine_mem_size = 4
        String docker_image
        Int preemptible_tries
    }
    String output_name = basename(input_bam, ".bam") + ".validation_report"
    Int disk_size = ceil(size(input_bam, "GB")) + addtional_disk_size
    Int command_mem_size = machine_mem_size - 1
    command {
        java -Xmx~{command_mem_size}G -jar /usr/gitc/picard.jar \
        ValidateSamFile \
        INPUT=~{input_bam} \
        OUTPUT=~{output_name} \
        MODE=SUMMARY \
        IS_BISULFITE_SEQUENCED=false
    }
    #Run time attributes:
    #Use a docker with the picard.jar. Set this up as a workspace attribute.
    #Read more about return codes here: https://github.com/broadinstitute/cromwell#continueonreturncode
    runtime {
        docker: docker_image
        memory: machine_mem_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible_tries
        continueOnReturnCode: [0,1]
    }
    #A text file is generated that will list errors or warnings that apply.
    output {
        File report = "~{output_name}"
    }
}

task AddOrReplaceReadGroups {
    input {
        File inputBam
        String sample_name
        String readgroupLibrary
        String readgroupPlatform
        String readgroupRunBarcode
        Int memoryGb
        Int diskSpaceGb
    }
    command <<<

        java -jar /usr/gitc/picard.jar AddOrReplaceReadGroups \
            I=~{inputBam} \
            O=~{sample_name}.readgroupadded.bam \
            RGID=1 \
            RGLB=~{readgroupLibrary} \
            RGPL=~{readgroupPlatform} \
            RGPU=~{readgroupRunBarcode} \
            RGSM=~{sample_name}

        samtools index ~{sample_name}.readgroupadded.bam
    >>>

    output {
        File bamWithReadGroupAdded = "${sample_name}.readgroupadded.bam"
        File bamWithReadGroupAddedIndex = "${sample_name}.readgroupadded.bam.bai"
    }

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
        memory: "${memoryGb} GB"
        cpu: "1"
        disks: "local-disk ${diskSpaceGb} HDD"
    }
}