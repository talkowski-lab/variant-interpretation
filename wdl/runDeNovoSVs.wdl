version 1.0
    
import "https://raw.githubusercontent.com/broadinstitute/gatk-sv/v0.26-beta/wdl/Structs.wdl"

workflow deNovoSV {

    input {

        File ped_input
        File python_config
        File vcf_file
        Array[String] contigs
        File genomic_disorder_input
        File raw_files_list
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override

    }

    call getGenomicDisorders{
        input:
            genomic_disorder_input=genomic_disorder_input,
            vcf_file=vcf_file,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_override
    }


    Array[String] raw_files = transpose(read_tsv(raw_files_list))[0]

    scatter(raw_file in raw_files){
        call raw_VcfToBed {
            input:
                vcf_file=raw_file,
                variant_interpretation_docker=variant_interpretation_docker,
                prefix = basename(raw_file, ".vcf.gz"),
                runtime_attr_override = runtime_attr_override
        }
    }

    call raw_mergeBed {
            input:
                bed_files=raw_VcfToBed.bed_output,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_override
    }

    scatter (contig in contigs){
        call subsetVcf {
            input:
                vcf_file=vcf_file,
                chromosome=contig,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_override
        }

        call vcfToBed{
            input:
                vcf_file=subsetVcf.vcf_output,
                chromosome=contig,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_override
        }    

        call raw_divideByChrom {
            input:
                bed_file = raw_mergeBed.concat_bed_output,
                chromosome = contig,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_override
        }

        call raw_reformatBed{
            input:
                per_chromosome_bed_file = raw_divideByChrom.per_chromosome_bed_output,
                ped_input=ped_input,
                chromosome=contig,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_override
        }

        call getDeNovo{
            input:
                bed_input=vcfToBed.bed_output,
                ped_input=ped_input,
                vcf_input=subsetVcf.no_header_vcf_output,
                disorder_input=getGenomicDisorders.gd_output,
                chromosome=contig,
                raw_proband=raw_reformatBed.reformatted_proband_output,
                raw_parents=raw_reformatBed.reformatted_parents_output,
                python_config=python_config,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_override
        }
    }

    call plot_mergeFinalBedFiles{
        input:
            bed_files = getDeNovo.denovo_output,
            outliers_files = getDeNovo.denovo_outliers,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_override
    }

    call plot_createPlots{
        input:
            bed_file = plot_mergeFinalBedFiles.concat_final_bed_output,
            outliers_file = plot_mergeFinalBedFiles.concat_final_bed_outliers_output,
            ped_input = ped_input,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_override
    }


    output {
    
        File concat_output = plot_mergeFinalBedFiles.concat_final_bed_output
        File concat_outlier_output = plot_mergeFinalBedFiles.concat_final_bed_outliers_output
        File denovo_output_plots = plot_createPlots.output_plots

    }
}

task getDeNovo{
    input{
        File bed_input
        File ped_input
        File vcf_input
        File disorder_input
        String chromosome
        File raw_proband
        File raw_parents
        File python_config
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File denovo_output = "~{chromosome}.denovo.bed"
        File denovo_outliers = "~{chromosome}.denovo.outliers.bed"
    }

    command <<<
            python3.9 /src/variant-interpretation/scripts/deNovoSVs.py \
                --bed ~{bed_input} \
                --ped ~{ped_input} \
                --vcf ~{vcf_input} \
                --disorder ~{disorder_input} \
                --out ~{chromosome}.denovo.bed \
                --raw_proband ~{raw_proband} \
                --raw_parents ~{raw_parents} \
                --config ~{python_config} \
                --verbose True \
                --outliers ~{chromosome}.denovo.outliers.bed
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}


task subsetVcf{
    input{
        File vcf_file
        String chromosome
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible_tries: 3,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File vcf_output = "~{chromosome}.vcf.gz"
        File no_header_vcf_output = "~{chromosome}.noheader.vcf.gz"
    }

    command <<<
        bcftools index ~{vcf_file}

        bcftools view ~{vcf_file} --regions ~{chromosome} -O z -o  ~{chromosome}.vcf.gz

        bcftools view ~{chromosome}.vcf.gz | grep -v ^## | bgzip -c > ~{chromosome}.noheader.vcf.gz

        
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task vcfToBed{
    input{
        File vcf_file
        String chromosome
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible_tries: 3,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File bed_output = "~{chromosome}.bed.gz"
    }

    command <<<
        svtk vcf2bed ~{vcf_file} --info ALL --include-filters ~{chromosome}.bed
        bgzip ~{chromosome}.bed
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task getGenomicDisorders{
    input{
        File vcf_file
        File genomic_disorder_input
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible_tries: 3,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File gd_output = "annotated.gd.variants.names.txt"
    }

    command <<<
        bedtools intersect -wa -wb -f 0.5 -a ~{vcf_file} -b ~{genomic_disorder_input} | cut -f 3 |sort -u> annotated.gd.variants.names.txt
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task raw_VcfToBed{
    input{
        File vcf_file
        String variant_interpretation_docker
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible_tries: 3,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File bed_output = "${prefix}.bed"
    }

    command {

        #convert from vcf to bed file
        
        svtk vcf2bed ~{vcf_file} --info SVTYPE ${prefix}.bed

    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
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

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible_tries: 3,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File concat_bed_output = "concat.bed"
    }

    command {

        cat ${sep=" " bed_files} > concat.bed

    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
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

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible_tries: 3,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File per_chromosome_bed_output = "${chromosome}.bed"
    }

    command {

         grep -w ^${chromosome} ${bed_file} > ${chromosome}.bed

    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
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

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible_tries: 3,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
       File reformatted_proband_output = "${chromosome}.proband.reformatted.sorted.bed"
       File reformatted_parents_output = "${chromosome}.parents.reformatted.sorted.bed"
    }

    command {

        #reformat bed file
        Rscript /src/variant-interpretation/scripts/reformatRawBed.R ${per_chromosome_bed_file} ${ped_input} ${chromosome}.proband.reformatted.bed ${chromosome}.parents.reformatted.bed
        sortBed -i ${chromosome}.proband.reformatted.bed > ${chromosome}.proband.reformatted.sorted.bed
        sortBed -i ${chromosome}.parents.reformatted.bed > ${chromosome}.parents.reformatted.sorted.bed

    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}   

task plot_mergeFinalBedFiles{
    input{
        Array[File] bed_files
        Array[File] outliers_files
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible_tries: 3,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File concat_final_bed_output = "final.denovo.merged.bed.gz"
        File concat_final_bed_outliers_output = "final.denovo.outliers.merged.bed.gz"
    }

    command {
        head -n+1 ${bed_files[1]} > final.denovo.merged.bed
        tail -n+2 -q ${sep=" " bed_files} >> final.denovo.merged.bed
        bgzip final.denovo.merged.bed
        
        head -n+1 ${outliers_files[1]} > final.denovo.outliers.merged.bed
        tail -n+2 -q ${sep=" " outliers_files} >> final.denovo.outliers.merged.bed
        bgzip final.denovov.outliers.merged.bed
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task plot_createPlots{
    input{
        File bed_file
        File outliers_file
        File ped_input
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible_tries: 3,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File output_plots = "output_plots.pdf"
    }

    command {

        Rscript /src/variant-interpretation/scripts/denovoSV_plots.R ${bed_file} ${outliers_file} ${ped_input} output_plots.pdf

    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}