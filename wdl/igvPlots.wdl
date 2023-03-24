version 1.0

##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

import "Structs.wdl"

workflow IGV {
    input{
        File varfile
        File Fasta
        File Fasta_idx
        File Fasta_dict
        File nested_repeats
        File simple_repeats
        File empty_track
        File ped_file
        Array[File] sample_crai_cram
        String buffer
        String buffer_large
        String igv_docker
        RuntimeAttr? runtime_attr_igv
    }

    scatter (file in sample_crai_cram) {
        Array[String] samples = transpose(read_tsv(file))[0]
        Array[File] crais = transpose(read_tsv(file))[1]
        Array[File] crams = transpose(read_tsv(file))[2]
        call runIGV_whole_genome{
            input:
                varfile = varfile,
                fasta = Fasta,
                fasta_dict = Fasta_dict,
                fasta_idx = Fasta_idx,
                nested_repeats = nested_repeats,
                simple_repeats = simple_repeats,
                empty_track = empty_track,
                ped_file = ped_file,
                samples = samples,
                crams = crams,
                crais = crais,
                sample_crai_cram = file,
                buffer = buffer,
                buffer_large = buffer_large,
                igv_docker = igv_docker,
                runtime_attr_override = runtime_attr_igv
        }
    }


    output{
        Array[File] tar_gz_pe = runIGV_whole_genome.pe_plots
    }
}

task runIGV_whole_genome{
    input{
        File varfile
        File fasta
        File fasta_idx
        File fasta_dict
        File nested_repeats
        File simple_repeats
        File empty_track
        File ped_file
        File sample_crai_cram
        Array[String] samples
        Array[File] crams
        Array[File] crais
        String buffer
        String buffer_large
        String igv_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([varfile, fasta, fasta_idx, fasta_dict, nested_repeats, simple_repeats, empty_track, ped_file, crams, crais]), "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: ceil(base_mem_gb),
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    File samples_file = write_lines(samples)
    File crams_file = write_lines(crams)
    String family = basename(sample_crai_cram, ".subset_sample_crai_cram.txt")
    command <<<
            set -euo pipefail
            #export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
            mkdir pe_igv_plots
            cat ~{varfile} | gunzip | cut -f1-6 > updated_varfile.bed
            grep -w -f ~{samples_file} updated_varfile.bed | cut -f1-5 | awk '{print $1,$2,$3,$4,$5}' | sed -e 's/ /\t/g' > ~{family}.bed
            python /src/makeigvpesr.py -v ~{family}.bed -n ~{nested_repeats} -s ~{simple_repeats} -e ~{empty_track} -f ~{fasta} -fam_id ~{family} -samples ~{samples_file} -crams ~{crams_file} -p ~{ped_file} -o pe_igv_plots -b ~{buffer} -l ~{buffer_large} -i pe.txt -bam pe.sh
            bash pe.sh
            xvfb-run --server-args="-screen 0, 1920x540x24" bash /IGV_2.4.14/igv.sh -b pe.txt
            tar -czf ~{family}_pe_igv_plots.tar.gz pe_igv_plots
        >>>
    
    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: igv_docker
    }
    output{
        File pe_plots= "~{family}_pe_igv_plots.tar.gz"
        File pe_txt = "pe.txt"
        File per_family_bed = "~{family}.bed"
        }
    }