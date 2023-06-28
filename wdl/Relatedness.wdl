version 1.0

import "Structs.wdl"

workflow relatedness {

    input {
        File vcf_list
        File ped_file
        File positions
        String relatedness_docker

        RuntimeAttr? runtime_attr_override_subset
        RuntimeAttr? runtime_attr_override_merge
        RuntimeAttr? runtime_attr_override_plink
    }

    Array[File] vcf_files = read_lines(vcf_list)

    scatter (vcf in vcf_files) {
        call subsetVCF{
            input:
                vcf_input=vcf,
                positions=positions,
                docker = relatedness_docker,
                runtime_attr_override = runtime_attr_override_subset
        }

        call mergeVCF{
            input:
                input_vcfs=subsetVCF.vcf_output,
                docker = relatedness_docker,
                runtime_attr_override = runtime_attr_override_merge
        }
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

task subsetVCF{
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
    }

    command <<<
        bcftools view -R ~{positions} ~{vcf_input} -O z -o ~{output_name}
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
        bcftools merge ~{sep=' ' input_vcfs} -Oz -o merged.5kpurcell.vcf.gz
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