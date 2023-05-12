version 1.0

##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

import "igvPlots.wdl" as igv_plots
import "Structs.wdl"

workflow IGV_all_samples {
    input {
        File ped_file
        File? fam_ids
        File sample_crai_cram
        File varfile
        File reference
        File reference_index
        Boolean cram_localization
        String prefix
        String buffer
        String buffer_large
        String sv_base_mini_docker
        String igv_docker
        RuntimeAttr? runtime_attr_run_igv
        RuntimeAttr? runtime_attr_igv
    }

    if (defined(fam_ids)) {
        File fam_ids_ = select_first([fam_ids])
        Array[String] family_ids = transpose(read_tsv(fam_ids_))[0]
    }

    if (!(defined(fam_ids))) {
        call generate_families{
            input:
                varfile = varfile,
                ped_file = ped_file,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_run_igv
        }
    }
    scatter (family in select_first([family_ids, generate_families.families])){
        call generate_per_family_sample_crai_cram{
            input:
                family = family,
                ped_file = ped_file,
                sample_crai_cram = sample_crai_cram,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_run_igv
            }
        call generate_per_family_bed{
            input:
                varfile = varfile,
                samples = generate_per_family_sample_crai_cram.per_family_samples,
                family = family,
                ped_file = ped_file,
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_run_igv
        }
        
        if (cram_localization){
            call igv_plots.IGV as IGV_localize {
                input:
                    varfile = generate_per_family_bed.per_family_varfile,
                    family = family,
                    ped_file = ped_file,
                    samples = generate_per_family_sample_crai_cram.per_family_samples,
                    crams_localize = generate_per_family_sample_crai_cram.per_family_crams_files,
                    crams_localize = generate_per_family_sample_crai_cram.per_family_crais_files,
                    buffer = buffer,
                    buffer_large = buffer_large,
                    reference = reference,
                    reference_index = reference_index,
                    igv_docker = igv_docker,
                    runtime_attr_igv = runtime_attr_igv
            }
        }

        if (!(cram_localization)){
            call igv_plots.IGV as IGV_parse {
                input:
                    varfile = generate_per_family_bed.per_family_varfile,
                    family = family,
                    ped_file = ped_file,
                    samples = generate_per_family_sample_crai_cram.per_family_samples,
                    crams_parse = generate_per_family_sample_crai_cram.per_family_crams_strings,
                    crams_parse = generate_per_family_sample_crai_cram.per_family_crais_strings,
                    buffer = buffer,
                    buffer_large = buffer_large,
                    reference = reference,
                    reference_index = reference_index,
                    igv_docker = igv_docker,
                    runtime_attr_igv = runtime_attr_igv
            }
        }
    }
    
    call integrate_igv_plots{
        input:
            igv_tar = select_first([IGV_localize.tar_gz_pe,IGV_parse.tar_gz_pe]),
            prefix = prefix, 
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_run_igv
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
    Float input_size = size(select_all([varfile, ped_file]), "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + input_size),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    command <<<
        set -euo pipefail
        cat ~{varfile} | gunzip | tail -n+2 | cut -f6 | tr ',' '\n' | sort -u > samples.txt #must have header line
        grep -w -f samples.txt ~{ped_file} | cut -f1 | sort -u  > families.txt
        >>>

    output{
        Array[String] families = read_lines("families.txt")
    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: sv_base_mini_docker
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
    Float input_size = size(select_all([sample_crai_cram, ped_file]), "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + input_size),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail
        grep -w ~{family} ~{ped_file} | cut -f2 > samples_list.txt
        grep -f samples_list.txt ~{sample_crai_cram} > subset_sample_crai_cram.txt
        cut -f1 subset_sample_crai_cram.txt > samples.txt
        cut -f2 subset_sample_crai_cram.txt > crai.txt
        cut -f3 subset_sample_crai_cram.txt > cram.txt
        >>>

    output{
        Array[String] per_family_samples = read_lines("samples.txt")
        Array[File] per_family_crams_files = read_lines("cram.txt")
        Array[File] per_family_crais_files = read_lines("crai.txt")
        Array[String] per_family_crams_strings = read_lines("cram.txt")
        Array[String] per_family_crais_strings = read_lines("crai.txt")
    }

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
    Float input_size = size(select_all([varfile, ped_file]), "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + input_size),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    String filename = basename(varfile, ".bed")

    command <<<
        set -euo pipefail
        cat ~{varfile} | gunzip | cut -f1-6 > updated_varfile.bed
        grep -f ~{write_lines(samples)} updated_varfile.bed | cut -f1-5 | awk '{print $1,$2,$3,$4,$5}' | sed -e 's/ /\t/g' > ~{filename}.~{family}.bed
        >>>

    output{
        File per_family_varfile = "~{filename}.~{family}.bed"
        }

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
    Float input_size = size(igv_tar, "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + input_size),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
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