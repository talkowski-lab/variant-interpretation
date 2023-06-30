#Author: Alba Sanchis-Juan <asanchis@broadinstitute.org>
#This script runs king relatedness on a VCF given a pedigree file

version 1.0

import "Structs.wdl"

workflow relatedness {

    input {
        File vcf_list
        File ped_file
        File positions
        Boolean split_by_chr
        String relatedness_docker

        Array[String] contigs

        RuntimeAttr? runtime_attr_override_subset
        RuntimeAttr? runtime_attr_override_merge
        RuntimeAttr? runtime_attr_override_plink
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

        call subsetPositionsVCF{
            input:
                vcf_input=vcf,
                positions=positions,
                docker = relatedness_docker,
                runtime_attr_override = runtime_attr_override_subset
        }

    }

    call mergeVCF{
        input:
            input_vcfs=subsetPositionsVCF.vcf_output,
            input_vcfs_index=subsetPositionsVCF.vcf_output_index,
            docker = relatedness_docker,
            runtime_attr_override = runtime_attr_override_merge
    }

    call runPlink{
        input:
            input_vcf = mergeVCF.merged_vcf,
            ped_file = ped_file,
            docker = relatedness_docker,
            runtime_attr_override = runtime_attr_override_plink
    }

    output{
        File kin = runPlink.kin
        File kin0 = runPlink.kin0
        File smiss = runPlink.smissing
        File vmiss = runPlink.vmissing
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

    String output_name = basename(vcf_input, "vcf.gz") + "5kpurcell.vcf.gz"

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
        File merged_vcf = "merged.5kpurcell.vcf.gz"
    }

    command <<<
        bcftools concat ~{sep=' ' input_vcfs} -Oz -o merged.5kpurcell.vcf.gz
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

task runPlink{
    input{
        File input_vcf
        File ped_file
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
        File smissing = "missing_output.plink.smiss"
        File vmissing = "missing_output.plink.vmiss"
        File kin = "king.kin"
        File kin0 = "king.kin0"

    }

    command <<<
        ##Prepare for plink
        plink --vcf ~{input_vcf} --maf 0.01 --make-bed --out output.plink

        ##All variants need an ID
        awk -F'\t' 'BEGIN {OFS = FS} $2=="." {print $1, $1 "_" $4, $3, $4, $5, $6} $2!="." {print $1,$2,$3,$4,$5,$6}' output.plink.bim > output.plink.bim.new
        mv output.plink.bim ori.output.plink.bim.new
        mv output.plink.bim.new output.plink.bim

        ##Check missing rate
        plink --bfile output.plink --missing --out missing_output.plink

        ##Change pedigree file to fam file
        mv ~{ped_file} output.plink.fam

        ##Run king
        king -b output.plink.bed --kinship --degree 2
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