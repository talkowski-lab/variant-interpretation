version 1.0

## Copyright Broad Institute, 2019
##
## Adaptation of multiple worflows to process CRAM to individual VCF
## Author: asanchis@broadinstitute.org
##
## Adapted from:
## CRAM to BAM: github.com/gatk-workflows/seq-format-conversion/CRAM-to-BAM:master
## addOrReplaceReadGroups: margolis/addOrReplaceReadGroups/1
## Haplotype Caller: github.com/gatk-workflows/gatk4-germline-snps-indels/haplotypecaller-gvcf-gatk4:2.3.1
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.

# WORKFLOW DEFINITION

workflow CramToHaplotypeCallerDragenFlow {
    input {
        File ref_fasta
        File ref_fasta_index_fai
        File ref_fasta_index_gzi
        File ref_dict
        File bam_or_cram
        String sample_name
        String gotc_docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
        Int preemptible_tries = 3
        File scattered_calling_intervals_list
        String readgroup_library
        String readgroup_platform
        String readgroup_run_barcode
        Int memory_addOrReplaceGroups
        Int disk_addOrReplaceGroups
        Boolean make_gvcf = true
        Boolean make_bamout = false
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.9.0"
        String gatk_path = "/gatk/gatk"
        String gitc_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
        String samtools_path = "samtools"
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

    # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
    # If we take the number we are scattering by and reduce by 20 we will have enough disk space
    # to account for the fact that the data is quite uneven across the shards.
    Array[File] scattered_calling_intervals = read_lines(scattered_calling_intervals_list)
    Int potential_hc_divisor = length(scattered_calling_intervals) - 20
    Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

    #Define output filename
    String sample_basename = basename(AddOrReplaceReadGroups.bamWithReadGroupAdded, ".bam")
    String vcf_basename = sample_basename
    String output_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
    String output_filename = vcf_basename + output_suffix

    # Call variants in parallel over grouped calling intervals
    scatter (interval_file in scattered_calling_intervals) {

        # Generate GVCF by interval
        call HaplotypeCaller {
            input:
                input_bam = AddOrReplaceReadGroups.bamWithReadGroupAdded,
                input_bam_index = AddOrReplaceReadGroups.bamWithReadGroupAddedIndex,
                interval_list = interval_file,
                output_filename = output_filename,
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index_fai = ref_fasta_index_fai,
                ref_fasta_index_gzi = ref_fasta_index_gzi,
                hc_scatter = hc_divisor,
                make_gvcf = make_gvcf,
                make_bamout = make_bamout,
                docker = gatk_docker,
                gatk_path = gatk_path
        }
    }

    # Merge per-interval GVCFs
    call MergeGVCFs {
        input:
            input_vcfs = HaplotypeCaller.output_vcf,
            input_vcfs_indexes = HaplotypeCaller.output_vcf_index,
            output_filename = output_filename,
            docker = gatk_docker,
            gatk_path = gatk_path
    }

    # Outputs that will be retained when execution is complete
    output {
        File output_vcf = MergeGVCFs.output_vcf
        File output_vcf_index = MergeGVCFs.output_vcf_index
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

        samtools view -h ~{inputBam} | \
            sed -e 's/DRAGEN_RGID/0/g' | \
            samtools view -h -O bam -o tmp.bam

        java -jar /usr/gitc/picard.jar AddOrReplaceReadGroups \
            I=tmp.bam \
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

# HaplotypeCaller per-sample in GVCF mode

task HaplotypeCaller {
    input {
        # Command parameters
        File input_bam
        File input_bam_index
        File interval_list
        String output_filename
        File ref_dict
        File ref_fasta
        File ref_fasta_index_fai
        File ref_fasta_index_gzi
        Float? contamination
        Boolean make_gvcf
        Boolean make_bamout
        Int hc_scatter

        String? gcs_project_for_requester_pays

        String gatk_path
        String? java_options

        # Runtime parameters
        String docker
        Int? mem_gb
        Int? disk_space_gb
        Boolean use_ssd = false
        Int? preemptible_attempts
    }

    String java_opt = select_first([java_options, "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"])

    Int machine_mem_gb = select_first([mem_gb, 7])
    Int command_mem_gb = machine_mem_gb - 1

    Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index_fai, "GB") + size(ref_dict, "GB")
    Int disk_size = ceil(((size(input_bam, "GB") + 30) / hc_scatter) + ref_size) + 20

    String vcf_basename = if make_gvcf then  basename(output_filename, ".gvcf") else basename(output_filename, ".vcf")
    String bamout_arg = if make_bamout then "-bamout ~{vcf_basename}.bamout.bam" else ""

    parameter_meta {
        input_bam: {
                       description: "a bam file",
                       localization_optional: false
                   }
        input_bam_index: {
                             description: "an index file for the bam input",
                             localization_optional: false
                         }
    }
    command {
        set -e

        ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G ~{java_opt}" \
        HaplotypeCaller \
        -R ~{ref_fasta} \
        -I ~{input_bam} \
        -L ~{interval_list} \
        -O ~{output_filename} \
        -contamination ~{default="0" contamination} \
        -G StandardAnnotation -G StandardHCAnnotation ~{true="-G AS_StandardAnnotation" false="" make_gvcf} \
        -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
        ~{true="-ERC GVCF" false="" make_gvcf} \
        ~{if defined(gcs_project_for_requester_pays) then "--gcs-project-for-requester-pays ~{gcs_project_for_requester_pays}" else ""} \
        ~{bamout_arg}

        # Cromwell doesn't like optional task outputs, so we have to touch this file.
        touch ~{vcf_basename}.bamout.bam
    }
    runtime {
        docker: docker
        memory: machine_mem_gb + " GB"
        disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 3])
    }
    output {
        File output_vcf = "~{output_filename}"
        File output_vcf_index = "~{output_filename}.tbi"
        File bamout = "~{vcf_basename}.bamout.bam"
    }
}
# Merge GVCFs generated per-interval for the same sample

task MergeGVCFs {
    input {
        # Command parameters
        Array[File] input_vcfs
        Array[File] input_vcfs_indexes
        String output_filename

        String gatk_path

        # Runtime parameters
        String docker
        Int? mem_gb
        Int? disk_space_gb
        Int? preemptible_attempts
    }
    Boolean use_ssd = false
    Int machine_mem_gb = select_first([mem_gb, 3])
    Int command_mem_gb = machine_mem_gb - 1

    command {
        set -e

        ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G"  \
        MergeVcfs \
        --INPUT ~{sep=' --INPUT ' input_vcfs} \
        --OUTPUT ~{output_filename}
    }
    runtime {
        docker: docker
        memory: machine_mem_gb + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 3])
    }
    output {
        File output_vcf = "~{output_filename}"
        File output_vcf_index = "~{output_filename}.tbi"
    }
}