version 1.0

##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

import "igv_trio_plots.wdl" as igv
import "Structs.wdl"

workflow IGV_all_samples {
    input {
        Array[String] pb_list
        Array[String] fa_list
        Array[String] mo_list
        Array[File] pb_cram_list
        Array[File] fa_cram_list
        Array[File] mo_cram_list
        File varfile
        File Fasta
        File Fasta_dict
        File Fasta_idx
        String prefix
        String sv_base_mini_docker
        String igv_docker
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    scatter (i in range(length(pb_list))){
        call generate_per_sample_bed{
            input:
                varfile = varfile,
                sample_id = pb_list[i],
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_override
        }

        call x_createPbCrais{
                input:
                    cram_file = pb_cram_list[i],
                    variant_interpretation_docker=variant_interpretation_docker,
                    runtime_attr_override = runtime_attr_override
            }

        call x_createMoCrais{
                input:
                    cram_file = mo_cram_list[i],
                    variant_interpretation_docker=variant_interpretation_docker,
                    runtime_attr_override = runtime_attr_override
            }

        call x_createFaCrais{
                input:
                    cram_file = fa_cram_list[i],
                    variant_interpretation_docker=variant_interpretation_docker,
                    runtime_attr_override = runtime_attr_override
            }


        call igv.IGV_trio as IGV_trio {
            input:
                varfile=generate_per_sample_bed.per_sample_varfile,
                Fasta = Fasta,
                Fasta_idx = Fasta_idx,
                Fasta_dict = Fasta_dict,
                pb=pb_list[i],
                fa=fa_list[i],
                mo=mo_list[i],
                pb_cram=pb_cram_list[i],
                fa_cram=fa_cram_list[i],
                mo_cram=mo_cram_list[i],
                pb_crai=x_createPbCrais.crai_file,
                fa_crai=x_createMoCrais.crai_file,
                mo_crai=x_createFaCrais.crai_file,
                igv_docker = igv_docker
                }
        }
        
    call integrate_igv_plots{
        input:
            igv_tar = IGV_trio.tar_gz_pe,
            prefix = prefix, 
            sv_base_mini_docker = sv_base_mini_docker
    }

    output{
        File tar_gz_pe = integrate_igv_plots.plot_tar
    }
    }


task generate_per_sample_bed{
    input {
        File varfile
        String sample_id
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }
    RuntimeAttr default_attr=object {
        cpu: 1,
        mem_gb: 1,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible: 1,
        max_retries: 1
    }

    String filename = basename(varfile, ".bed")
    command <<<
        grep -w ~{sample_id} ~{varfile} | cut -f1-5 | awk '{print $1,$2,$3,$5,$4}' | sed -e 's/ /\t/g' > ~{filename}.~{sample_id}.bed
        >>>

    output{
        File per_sample_varfile= "~{filename}.~{sample_id}.bed"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_base_mini_docker
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    }

task integrate_igv_plots{
    input {
        Array[File] igv_tar
        String prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }
    RuntimeAttr default_attr=object {
        cpu: 1,
        mem_gb: 1,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible: 1,
        max_retries: 1
    }

    command <<<
        mkdir ~{prefix}_igv_plots
        while read file; do
            tar -zxf ${file}
            mv pe_igv_plots/*  ~{prefix}_igv_plots/
        done < ~{write_lines(igv_tar)};
        tar -czf ~{prefix}_igv_plots.tar.gz ~{prefix}_igv_plots
    >>>

    output{
        File plot_tar = "~{prefix}_igv_plots.tar.gz"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_base_mini_docker
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task x_createPbCrais{
    input{
        File cram_file
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

        File crai_file = "${cram_file}.crai"
    }


    command {

        samtools index -b ${cram_file}
    }

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

task x_createMoCrais{
    input{
        File cram_file
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

        File crai_file = "${cram_file}.crai"
    }


    command {

        samtools index -b ${cram_file}
    }

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

task x_createFaCrais{
    input{
        File cram_file
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

        File crai_file = "${cram_file}.crai"
    }


    command {

        samtools index -b ${cram_file}
    }

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


