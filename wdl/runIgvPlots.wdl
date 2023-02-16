version 1.0

##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

import "igvPlots.wdl" as igv
import "Structs.wdl"

workflow IGV_all_samples {
    input {
        File ped_file
        File sample_crai_cram
        File varfile
        File Fasta
        File Fasta_dict
        File Fasta_idx
        File nested_repeats
        File simple_repeats
        File empty_track
        String prefix
        String sv_base_mini_docker
        String igv_docker
        RuntimeAttr? runtime_attr_override
    }

    call generate_families{
        input:
        varfile = varfile,
        ped_file = ped_file,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_override
    }
    scatter (family in generate_families.families){
        call generate_per_family_sample_crai_cram{
            input:
                family = family,
                ped_file = ped_file,
                sample_cram = sample_cram,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_override
            }
        call generate_per_family_bed{
            input:
                varfile = varfile,
                samples = samples,
                family = family,
                ped_file = ped_file,
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_override
            }
        call igv.IGV_trio as IGV_trio {
            input:
                varfile=generate_per_family_bed.per_family_varfile,
                Fasta = Fasta,
                Fasta_idx = Fasta_idx,
                Fasta_dict = Fasta_dict,
                nested_repeats = nested_repeats,
                simple_repeats = simple_repeats,
                empty_track = empty_track,
                sample = sample,
                ped_file = ped_file,
                samples = generate_per_family_sample_crai_cram.per_family_samples,
                crams = generate_per_family_sample_crai_cram.per_family_crams,
                crais = generate_per_family_sample_crai_cram.per_family_crais,
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

task generate_families{
    input {
        File varfile
        File ped_file
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
        set -euo pipefail
        cat ~{varfile} | gunzip | tail -n+2 | cut -f6 | tr ',' '\n' | sort -u > samples.txt #must have header line
        grep -w -f samples.txt ~{ped_file} | cut -f1 | sort -u  > families.txt
        >>>

    output{
        Array[String] families = read_lines("families.txt")
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

task generate_per_family_sample_crai_cram{
    input {
        String family
        File ped_file
        File sample_crai_cram
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
        set -euo pipefail
        grep -w ~{family} ~{ped_file} | cut -f2 > samples.txt
        grep -f samples.txt ~{sample_crai_cram} > subset_sample_crai_cram.txt
        cut -f2 subset_sample_crai_cram.txt > crai.txt
        cut -f3 subset_sample_crai_cram.txt > cram.txt
        >>>

    output{
        Array[String] per_family_samples = read_lines("samples.txt")
        Array[File] per_family_crams = read_lines("cram.txt")
        Array[File] per_family_crais = read_lines("crai.txt")
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

task generate_per_family_bed{
    input {
        File varfile
        Array[String] samples
        String family
        File ped_file
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
        set -euo pipefail
        cat ~{varfile} | gunzip | cut -f1-6 > updated_varfile.bed
        grep -f ~{write_lines(samples)} updated_varfile.bed | cut -f1-5 | awk '{print $1,$2,$3,$4,$5}' | sed -e 's/ /\t/g' > ~{filename}.~{family}.bed
        >>>

    output{
        File per_family_varfile= "~{filename}.~{family}.bed"
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