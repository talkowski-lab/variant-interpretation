version 1.0
    
import "Structs.wdl"

workflow reformatRawFiles {

    input {
        Array[String] contigs
        File raw_files_list
        File ped_input
        String variant_interpretation_docker
        Boolean depth
        RuntimeAttr? runtime_attr_vcf_to_bed
        RuntimeAttr? runtime_attr_merge_bed
        RuntimeAttr? runtime_attr_divide_by_chrom
        RuntimeAttr? runtime_attr_reformat_bed
    }

    Array[String] raw_files = transpose(read_tsv(raw_files_list))[0]

    scatter(raw_file in raw_files){
        call raw_VcfToBed {
            input:
                vcf_file=raw_file,
                variant_interpretation_docker=variant_interpretation_docker,
                prefix = basename(raw_file, ".vcf.gz"),
                runtime_attr_override = runtime_attr_vcf_to_bed
        }
    }

    call raw_mergeBed {
            input:
                bed_files=raw_VcfToBed.bed_output,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_merge_bed
    }

    scatter (contig in contigs){
        call raw_divideByChrom {
            input:
                bed_file = raw_mergeBed.concat_bed_output,
                chromosome = contig,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_divide_by_chrom
        }

        if (depth) {
            call raw_reformatBedDepth{
                input:
                    per_chromosome_bed_file = raw_divideByChrom.per_chromosome_bed_output,
                    ped_input=ped_input,
                    chromosome=contig,
                    variant_interpretation_docker=variant_interpretation_docker,
                    runtime_attr_override = runtime_attr_reformat_bed
            }
        }

        else {
            call raw_reformatBed{
                input:
                    per_chromosome_bed_file = raw_divideByChrom.per_chromosome_bed_output,
                    ped_input=ped_input,
                    chromosome=contig,
                    variant_interpretation_docker=variant_interpretation_docker,
                    runtime_attr_override = runtime_attr_reformat_bed
            }
        }
    }


    output {
        Array[File] reformatted_parents_raw_files = select_first([raw_reformatBed.reformatted_parents_output, raw_reformatBedDepth.reformatted_parents_output])
        Array[File] reformatted_proband_raw_files = select_first([raw_reformatBed.reformatted_proband_output, raw_reformatBedDepth.reformatted_proband_output])
    }
}   

task raw_VcfToBed{
    input{
        File vcf_file
        String variant_interpretation_docker
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb + input_size * 3.0,
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File bed_output = "${prefix}.bed.gz"
    }

    command {
        set -euo pipefail

        #convert from vcf to bed file
        svtk vcf2bed ~{vcf_file} --info SVTYPE ${prefix}.bed
        bgzip ${prefix}.bed

    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}


task raw_mergeBed{
    input{
        Array[File] bed_files
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(bed_files, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb + input_size * 3.0,
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File concat_bed_output = "concat.bed.gz"
    }

    command {
        set -euo pipefail
        zcat ${sep=" " bed_files} | bgzip -c > concat.bed.gz

    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task raw_divideByChrom{
    input{
        File bed_file
        String chromosome
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(bed_file, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb + input_size * 3.0,
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File per_chromosome_bed_output = "${chromosome}.bed.gz"
    }

    command {
        set -euo pipefail
        
        zcat ${bed_file} | \
        grep -w ^${chromosome} | \
        bgzip -c > ${chromosome}.bed.gz

    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task raw_reformatBed{
    input{
        File per_chromosome_bed_file
        File ped_input
        String chromosome
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([per_chromosome_bed_file, ped_input]), "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb + input_size * 3.0,
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
       File reformatted_proband_output = "${chromosome}.proband.reformatted.sorted.bed.gz"
       File reformatted_parents_output = "${chromosome}.parents.reformatted.sorted.bed.gz"
    }

    command {
        set -euo pipefail

        #reformat bed file
        Rscript /src/variant-interpretation/scripts/reformatRawBed.R ${per_chromosome_bed_file} ${ped_input} ${chromosome}.proband.reformatted.bed ${chromosome}.parents.reformatted.bed
        sortBed -i ${chromosome}.proband.reformatted.bed | bgzip -c > ${chromosome}.proband.reformatted.sorted.bed.gz
        sortBed -i ${chromosome}.parents.reformatted.bed | bgzip -c > ${chromosome}.parents.reformatted.sorted.bed.gz

    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}   

task raw_reformatBedDepth{
    input{
        File per_chromosome_bed_file
        File ped_input
        String chromosome
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([per_chromosome_bed_file, ped_input]), "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb + input_size * 3.0,
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
       File reformatted_proband_output = "${chromosome}.proband.depth.reformatted.sorted.bed.gz"
       File reformatted_parents_output = "${chromosome}.parents.depth.reformatted.sorted.bed.gz"
    }

    command {
        set -euo pipefail

        #reformat bed file
        Rscript /src/variant-interpretation/scripts/reformatRawBed.R ${per_chromosome_bed_file} ${ped_input} ${chromosome}.proband.reformatted.bed ${chromosome}.parents.reformatted.bed
        sortBed -i ${chromosome}.proband.reformatted.bed | bgzip -c > ${chromosome}.proband.depth.reformatted.sorted.bed.gz
        sortBed -i ${chromosome}.parents.reformatted.bed | bgzip -c > ${chromosome}.parents.depth.reformatted.sorted.bed.gz

    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}   



