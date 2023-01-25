version 1.0
    
import "Structs.wdl"

workflow familyFiltering {

    input {
        File vcf_file
        File ped_file
        File genomic_disorder_input
        String cohort_prefix
        Array[String] families
        File rconfig
        File rconfig
        String variant_interpretation_docker

        RuntimeAttr? runtime_attr_override_vcfToBed
        RuntimeAttr? runtime_attr_override_getGD
        RuntimeAttr? runtime_attr_override_subsetFamily
        RuntimeAttr? runtime_attr_override_svFiltering
    }

    call vcfToBed{
        input:
            vcf_file=vcf_file,
            cohort_prefix = cohort_prefix,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_override_vcfToBed
    }
    call getGenomicDisorders{
        input:
            genomic_disorder_input=genomic_disorder_input,
            bed_file=vcfToBed.bed_output,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_override_getGD
    }

    scatter (family in families) {
        call subsetFamilyVCF{
            input:
                family=family,
                vcf_file=vcf_file,
                ped_file=ped_file,
                genomic_disorder_input=genomic_disorder_input,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_override_subsetFamily
        }

        call SVfamilyFiltering{
            input:
                family=family,
                family_vcf=subsetFamilyVCF.vcf_family,
                bed_file=vcfToBed.bed_output,
                ped_file=ped_file,
                genomic_disorder_names=getGenomicDisorders.gd_output,
                rconfig=rconfig,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_override_svFiltering
        }
    }
}

task vcfToBed{
    input{
        File vcf_file
        String cohort_prefix
        String variant_interpretation_docker
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
        File bed_output = "~{cohort_prefix}.bed.gz"
    }

    command <<<
        svtk vcf2bed ~{vcf_file} --info ALL --include-filters ~{cohort_prefix}.bed.gz
        bgzip ~{cohort_prefix}.bed.gz
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task getGenomicDisorders{
    input{
        File bed_file
        File genomic_disorder_input
        String variant_interpretation_docker
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
        File gd_output = "gd.names.txt"
    }

    command <<<
        bedtools intersect -wa -wb -f 0.5 -a ~{bed_file} -b ~{genomic_disorder_input} | cut -f 3 |sort -u> gd.names.txt
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task subsetFamilyVCF{
    input{
        String family
        File vcf_file
        File ped_file
        File genomic_disorder_input
        String variant_interpretation_docker
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
        File vcf_family = "~{family}.vcf.gz"
    }

    command <<<
        samples=`grep -w ^~{family} ~{ped_file} | cut -f2| sort -u |xargs| sed -e 's/ /,/g'`
        bcftools view ~{vcf_file} -s $samples --force-samples -O z -o ~{family}.vcf.gz
        bcftools index ~{family}.vcf.gz
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task SVfamilyFiltering{
    input{
        String family
        File family_vcf
        File bed_file
        File ped_file
        File genomic_disorder_names
        File rconfig
        String variant_interpretation_docker
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
        File family_filtered_svs = "~{family}.txt.gz"
    }

    command <<<
        Rscript familyFiltering.R \
            -f ~{family} \
            -g ~{family_vcf} \
            -i ~{bed_file} \
            -m ~{ped_file} \
            -d ~{genomic_disorder_names} \
            -c ~{rconfig} \
            -u /src/variant-interpretation/scripts/familyFilteringFunctions.R \
            -v
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}