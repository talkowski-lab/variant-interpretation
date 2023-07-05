#Author: Alba Sanchis-Juan <asanchis@broadinstitute.org>
#This script runs king relatedness on a VCF given a pedigree file

version 1.0

import "Structs.wdl"

workflow relatedness {

    input {
        File vcf_list
        File positions
        File exclude_input

        Boolean split_by_chr
        Boolean is_snv_indel

        String relatedness_docker

        Array[String] contigs

        RuntimeAttr? runtime_attr_override_subset
        RuntimeAttr? runtime_attr_override_merge
        RuntimeAttr? runtime_attr_override_vcftools
        RuntimeAttr? runtime_attr_override_split_by_chr
    }

    ##Split chromosome
    if (split_by_chr){

        File vcf_file = read_lines(vcf_list)[0]
        File vcf_file_index = vcf_file + ".tbi"

        scatter (contig in contigs) {

            call splitVCF{
                input:
                    input_vcf = vcf_file,
                    input_vcf_index = vcf_file_index,
                    contig = contig,
                    docker = relatedness_docker,
                    runtime_attr_override = runtime_attr_override_split_by_chr
            }
        }
    }

    ##If split by chr is set to false, just use the list of sharded VCFs
    Array[File] vcf_files = select_first([splitVCF.vcf_output, read_lines(vcf_list)])

    scatter (vcf in vcf_files) {

        File vcf_index = vcf + ".tbi"

        if (is_snv_indel) {
            call subsetPositionsVCF{
                input:
                    vcf_input=vcf,
                    positions=positions,
                    docker = relatedness_docker,
                    runtime_attr_override = runtime_attr_override_subset
            }
        }

        if (!is_snv_indel) {
            call subsetSVs{
                input:
                    vcf_input=vcf,
                    docker = relatedness_docker,
                    runtime_attr_override = runtime_attr_override_subset
            }

            call excludeSVregions{
                input:
                    vcf_input=vcf,
                    exclude_input=exclude_input,
                    docker = relatedness_docker,
                    runtime_attr_override = runtime_attr_override_subset
                }
            }
        }
    }

    Array[File] subset_files = select_first([subsetPositionsVCF.vcf_output, excludeSVregions.vcf_output, subsetSVs.vcf_output])
    Array[File] subset_files_index = select_first([subsetPositionsVCF.vcf_output_index, excludeSVregions.vcf_output_index, subsetSVs.vcf_output_index])

    call mergeVCF{
        input:
            input_vcfs=subset_files,
            input_vcfs_index=subset_files_index,
            docker = relatedness_docker,
            runtime_attr_override = runtime_attr_override_merge
    }

    call runVCFTools{
        input:
            input_vcf = mergeVCF.merged_vcf,
            docker = relatedness_docker,
            runtime_attr_override = runtime_attr_override_vcftools
    }

    output{
        File output_relatedness = runVCFTools.out_relatedness
    }
}


task splitVCF{
    input{
        File input_vcf
        File input_vcf_index
        String contig
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible: 3,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String output_name = basename(input_vcf, "vcf.gz") + contig + ".vcf.gz"

    output{
        File vcf_output = output_name
    }

    command <<<
        tabix -p vcf ~{input_vcf}
        bcftools view -r ~{contig} ~{input_vcf} -O z -o ~{output_name}
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}

task subsetPositionsVCF{
    input{
        File vcf_input
        File positions
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible: 3,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String output_name = basename(vcf_input, "vcf.gz") + "subset.vcf.gz"

    output{
        File vcf_output = output_name
        File vcf_output_index = output_name + ".tbi"
    }

    command <<<
        tabix -p vcf ~{vcf_input}
        bcftools view -R ~{positions} ~{vcf_input} -O z -o ~{output_name}
        tabix -p vcf ~{output_name}
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}

task subsetSVs{
    input{
        File vcf_input
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible: 3,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String output_name = basename(vcf_input, "vcf.gz") + "subset.vcf.gz"

    output{
        File vcf_output = output_name
        File vcf_output_index = output_name + ".tbi"
    }

    command <<<
        tabix -p vcf ~{vcf_input}
        bcftools view ~{vcf_input} | grep -E "^#|DEL|DUP|INS" | \
            bcftools view -i 'INFO/gnomad_v2.1_sv_AF>0.05' | \
            bcftools view -i 'AF>0.05' | \
            bcftools view -e 'AF>0.9' | \
            grep -vE "^chrX|^chrY" | \
            grep -v "HIGH_SR_BACKGROUND" | \
            grep -v "ALGORITHMS=wham;" | bgzip -c > ~{output_name}
        tabix -p vcf ~{output_name}
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}

task excludeSVregions{
    input{
        File vcf_input
        File exclude_input
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible: 3,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String output_name = basename(vcf_input, "vcf.gz") + "exclude.vcf.gz"

    output{
        File vcf_output = output_name
        File vcf_output_index = output_name + ".tbi"
    }

    command <<<
        bcftools view -h ~{vcf_input} > tmp.vcf
        bedtools intersect -a ~{vcf_input} -b ~{exclude_input} -v -f 0.5 -r | bgzip -c >> tmp.vcf
        mv tmp.vcf ~{output_name}
        tabix -p vcf ~{output_name}
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}

task mergeVCF{
    input{
        Array [File] input_vcfs
        Array [File] input_vcfs_index
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible: 3,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File merged_vcf = "merged.subset.vcf.gz"
    }
    command <<<
        bcftools concat ~{sep=' ' input_vcfs} -Oz -o merged.subset.vcf.gz
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}

task runVCFTools{
    input{
        File input_vcf
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu: 1,
        mem_gb: 8,
        disk_gb: 4,
        boot_disk_gb: 4,
        preemptible: 3,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File out_relatedness = "out.relatedness2"
    }

    command <<<
        ##Run VCFtools
        vcftools --gzvcf ~{input_vcf} --relatedness2
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}