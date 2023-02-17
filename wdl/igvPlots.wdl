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
        String family
        File ped_file
        Array[String] samples
        Array[File] crams
        Array[File] crais
        String igv_docker
    }

    call runIGV_whole_genome{
        input:
            varfile = varfile,
            fasta = Fasta,
            fasta_dict = Fasta_dict,
            fasta_idx = Fasta_idx,
            nested_repeats = nested_repeats,
            simple_repeats = simple_repeats,
            empty_track = empty_track,
            family = family,
            ped_file = ped_file,
            samples = samples,
            crams = crams,
            crais = crais,
            igv_docker = igv_docker
    }

    output{
        File tar_gz_pe = runIGV_whole_genome.pe_plots
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
        String family
        File ped_file
        Array[String] samples
        Array[File] crams
        Array[File] crais
        String igv_docker
    }


    command <<<
            set -euo pipefail
            #export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
            mkdir pe_igv_plots
            python /src/makeigvpesr.py -v ~{varfile} -n ~{nested_repeats} -s ~{simple_repeats} -e ~{empty_track} -f ~{fasta} -fam_id ~{family} -samples ~{sep="," samples} -crams ~{sep="," crams} -p ~{ped_file} -o pe_igv_plots -b 500
            bash pe.sh
            xvfb-run --server-args="-screen 0, 1920x540x24" bash /IGV_2.4.14/igv.sh -b pe.txt
            tar -czf ~{family}_pe_igv_plots.tar.gz pe_igv_plots

        >>>
    runtime {
        docker: igv_docker
        preemptible: 1
        memory: "15 GB"
        disks: "local-disk 100 HDD"
        }
    output{
        File pe_plots="~{family}_pe_igv_plots.tar.gz"
        File pe_txt = "pe.txt"
        }
    }