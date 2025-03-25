version 1.0
    
import "Structs2.wdl"

workflow familyFiltering {

    input {
#        File bed_file
        File vcf_file
        File ped_file
        File genomic_disorder_input
        String cohort_prefix
        Boolean run_compound_het
        Array[String] families
        Array[String] AF_columns

        String variant_interpretation_docker

        File genelist
        File eo_file
        File prec_file
        File pli_file
        File hpo_db
        File mim_file

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
            bed_file=vcfToBed.bed_output_ref,
#            bed_file=bed_file,
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
                run_compound_het=run_compound_het,
#                bed_file=bed_file,
                ped_file=ped_file,
                genomic_disorder_names=getGenomicDisorders.gd_output,
                genelist=genelist,
                eo_file=eo_file,
                prec_file=prec_file,
                pli_file=pli_file,
                hpo_db=hpo_db,
                AF_columns=AF_columns,
                mim_file=mim_file,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override=runtime_attr_override_svFiltering
        }
    }
    output{
        Array[File] family_filtered_svs = SVfamilyFiltering.family_filtered_svs
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
        File bed_output_ref = "~{cohort_prefix}.ref.bed.gz"
    }

    command <<<
        set -euo pipefail

        echo "Starting svtk"
        svtk vcf2bed -i ALL --include-filters ~{vcf_file} - | bgzip -c > ~{cohort_prefix}.bed.gz
        echo "svtk finished"

        echo "Starting reformat of bed file"
        zcat ~{cohort_prefix}.bed.gz | \
            grep -E "DEL|DUP" | \
            awk '{print $1"_"$5"\t"$2"\t"$3"\t"$4"\t"$5}' | \
#            grep -v ^# | \
            bgzip -c > ~{cohort_prefix}.ref.bed.gz
        echo "Reformat finished"

#        svtk vcf2bed ~{vcf_file} --info ALL --include-filters ~{cohort_prefix}.bed.gz
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
        bedtools intersect -wa -f 0.3 -r -a ~{bed_file} -b ~{genomic_disorder_input} | cut -f 4 |sort -u> gd.names.txt
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
        Boolean run_compound_het
        File genelist
        File eo_file
        File prec_file
        File pli_file
        File hpo_db
        File mim_file

        Array[String] AF_columns

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
        File family_filtered_svs = "~{family}.filt.txt"
    }

    command <<<

        echo "
            genelist_path <- '~{genelist}'
            eo_path <- '~{eo_file}'
            prec_path <- '~{prec_file}'
            pli_path <- '~{pli_file}'
            hpodb_path <- '~{hpo_db}'
            mim_path <- '~{mim_file}'
            af_columns <- '~{sep="," AF_columns}'
            comp_het <- '~{run_compound_het}'
        " > config.R

        Rscript /scripts/variant-interpretation/scripts/familyFiltering.R \
            -f ~{family} \
            -g ~{family_vcf} \
            -i ~{bed_file} \
            -m ~{ped_file} \
            -d ~{genomic_disorder_names} \
            -c config.R \
            -u /scripts/variant-interpretation/scripts/familyFilteringFunctions.R \
            -o ~{family}.filt.txt \
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