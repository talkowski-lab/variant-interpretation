version 1.0

##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

import "igvPlots.wdl" as igv
import "Structs.wdl"
import "GetShardInputs.wdl" as GetShardInputs

workflow IGV_all_samples {
    input {
        File ped_file
        File? fam_ids
        File sample_crai_cram
        File varfile
        File Fasta
        File Fasta_dict
        File Fasta_idx
        File nested_repeats
        File simple_repeats
        File empty_track
        String prefix
        String buffer
        String buffer_large
        Int families_per_shard
        String sv_base_mini_docker
        String igv_docker
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_run_igv
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

    Int num_families = length(select_first([family_ids,generate_families.families]))
    Float num_families_float = num_families
    Int num_shards = ceil(num_families_float / families_per_shard)

    scatter (i in range(num_shards)) {
        call GetShardInputs.GetShardInputs as GetShardFamilies {
            input:
                items_per_shard = families_per_shard,
                shard_number = i,
                num_items = num_families,
                all_items = select_first([family_ids,generate_families.families])
    }
        call generate_per_family_sample_crai_cram{
            input:
                families = GetShardFamilies.shard_items,
                ped_file = ped_file,
                sample_crai_cram = sample_crai_cram,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_run_igv
            }

        call generate_per_family_bed{
            input:
                varfile = varfile,
                ped_file = ped_file,
                families = GetShardFamilies.shard_items,
                ped_file = ped_file,
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_run_igv
            }
        
        call igv.IGV as IGV {
            input:
                varfiles = generate_per_family_bed.per_family_varfile,
                Fasta = Fasta,
                Fasta_idx = Fasta_idx,
                Fasta_dict = Fasta_dict,
                nested_repeats = nested_repeats,
                simple_repeats = simple_repeats,
                empty_track = empty_track,
                sample_crai_cram=sample_crai_cram,
                families = GetShardFamilies.shard_items,
                ped_file = ped_file,
                samples = generate_per_family_sample_crai_cram.per_family_samples,
                crams = generate_per_family_sample_crai_cram.per_family_crams,
                crais = generate_per_family_sample_crai_cram.per_family_crais,
                buffer = buffer,
                buffer_large = buffer_large,
                igv_docker = igv_docker,
                runtime_attr_run_igv = runtime_attr_run_igv
        }
    }
    call integrate_igv_plots{
        input:
            igv_tar = flatten(IGV.tar_gz_pe),
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
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: ceil(base_mem_gb + input_size * 3.0),
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
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
        File families = "families.txt"
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

task clusterPed{
    input {
        File families
        File ped_file
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(select_all([families, ped_file]), "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: ceil(base_mem_gb + input_size * 3.0),
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    command <<<
        set -euo pipefail
        grep -w -f ~{families} ~{ped_file} > updated_ped.txt
        Rscript src/variant-interpretation/scripts/clusterPed.R updated_ped.txt
        cat subset_ped.txt | cut -f7 | sort -u > clusters.txt
        i=0
        while read -r line;
        do
            let "i=$i+1"
            grep -w "$line" clustered_ped.txt | cut -f1 | sort -u > family_ids.$i.txt 
        done<clusters.txt

        >>>

    output{
        File clusters = "clusters.txt"
        Array[Array[String]] family_clusters = read_lines("family_ids.$i.txt")
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


task generate_per_family_sample_crai_cram{
    input {
        Array[String] families
        File ped_file
        File sample_crai_cram
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(select_all([sample_crai_cram, ped_file]), "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: ceil(base_mem_gb + input_size * 3.0),
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail
        for family in ~{sep=' ' families}
        do
            grep -w "${family}" ~{ped_file} | cut -f2 > samples_list.$family.txt
            grep -f samples_list."${family}".txt ~{sample_crai_cram} > subset_sample_crai_cram.$family.txt
            cut -f1 subset_sample_crai_cram.$family.txt > samples.$family.txt
            cut -f2 subset_sample_crai_cram.$family.txt > crai.$family.txt
            cut -f3 subset_sample_crai_cram.$family.txt > cram.$family.txt
        done;
        >>>

    output{
        Array[Array[String]] per_family_samples = read_lines("samples.$family.txt")
        Array[Array[File]] per_family_crams = read_lines("cram.$family.txt")
        Array[Array[File]] per_family_crais = read_lines("crai.$family.txt")
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
        File ped_file
        Array[String] families
        File ped_file
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(select_all([varfile, ped_file]), "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: ceil(base_mem_gb + input_size * 3.0),
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
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
        for family in ~{sep=' ' families}
        do
            grep -w "${family}" ~{ped_file} | cut -f2 | sort -u > samples.$family.txt
            grep -w -f samples.$family.txt updated_varfile.bed | cut -f1-5 | awk '{print $1,$2,$3,$4,$5}' | sed -e 's/ /\t/g' > ~{filename}.$family.bed
        done;
        >>>

    output{
        Array[File] per_family_varfile = "~{filename}.$family.bed"
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
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: ceil(base_mem_gb + input_size * 3.0),
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
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