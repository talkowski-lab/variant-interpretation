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
        String? buffer
        String? buffer_large
        String igv_docker
        RuntimeAttr? runtime_attr_run_igv
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
            buffer = buffer,
            igv_docker = igv_docker,
            runtime_attr_override = runtime_attr_run_igv
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
        String? buffer
        String? buffer_large
        String igv_docker
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

    String buff = select_first([buffer, 500])
    String large_buff = select_first([buffer_large, 1000])
    command <<<
            set -euo pipefail
            #export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
            mkdir pe_igv_plots
        python /src/makeigvpesr.py -v ~{varfile} -n ~{nested_repeats} -s ~{simple_repeats} -e ~{empty_track} -f ~{fasta} -fam_id ~{family} -samples ~{sep="," samples} -crams ~{sep="," crams} -p ~{ped_file} -o pe_igv_plots -b ~{buff} -l ~{large_buff}
            bash pe.sh
            xvfb-run --server-args="-screen 0, 1920x540x24" bash /IGV_2.4.14/igv.sh -b pe.txt
            tar -czf ~{family}_pe_igv_plots.tar.gz pe_igv_plots

        >>>
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
    output{
        File pe_plots="~{family}_pe_igv_plots.tar.gz"
        File pe_txt = "pe.txt"
        }
    }