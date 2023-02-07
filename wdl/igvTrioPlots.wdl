version 1.0

##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

import "Structs.wdl"

workflow IGV_trio {
    input{
        File varfile
        File Fasta
        File Fasta_idx
        File Fasta_dict
        File nested_repeats
        File simple_repeats
        File empty_track
        String fam_id
        File ped_file
        File sample_cram
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
            fam_id = fam_id,
            ped_file = ped_file,
            sample_cram = sample_cram,
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
        String fam_id
        File ped_file
        File sample_cram
        String igv_docker
    }
    command <<<
            set -euo pipefail
            #export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

            grep ~{fam_id} ~{ped_file} > sample_ids.txt
            grep -f sample_ids.txt ~{sample_cram} > subset_sample_cram.txt
            python /src/makeigvpesr_trio.py ~{varfile} ~{nested_repeats} ~{simple_repeats} ~{empty_track} ~{fasta} ~{fam_id} subset_sample_cram.txt pe_igv_plots -b 500
            xvfb-run --server-args="-screen 0, 1920x540x24" bash /IGV_2.4.14/igv.sh -b pe.txt
            tar -czf ~{fam_id}_pe_igv_plots.tar.gz pe_igv_plots

        >>>
    runtime {
        docker: igv_docker
        preemptible: 1
        memory: "15 GB"
        disks: "local-disk 100 HDD"
        }
    output{
        File pe_plots="~{fam_id}_pe_igv_plots.tar.gz"
        File pe_txt = "pe.txt"
        }
    }