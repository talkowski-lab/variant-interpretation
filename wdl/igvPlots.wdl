version 1.0

##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

import "Structs.wdl"

workflow IGV {
    input{
        Array[File] varfiles
        File Fasta
        File Fasta_idx
        File Fasta_dict
        File nested_repeats
        File simple_repeats
        File empty_track
        File sample_crai_cram
        Array[String] families
        File ped_file
        Array[Array[String]] samples
        Array[Array[File]] crams
        Array[Array[File]] crais
        String buffer
        String buffer_large
        String igv_docker
        RuntimeAttr? runtime_attr_run_igv
    }

    call runIGV_whole_genome{
        input:
            varfiles = varfiles,
            fasta = Fasta,
            fasta_dict = Fasta_dict,
            fasta_idx = Fasta_idx,
            nested_repeats = nested_repeats,
            simple_repeats = simple_repeats,
            empty_track = empty_track,
            sample_crai_cram=sample_crai_cram,
            families = families,
            ped_file = ped_file,
            samples = samples,
            crams = crams,
            crais = crais,
            buffer = buffer,
            buffer_large = buffer_large,
            igv_docker = igv_docker,
            runtime_attr_override = runtime_attr_run_igv
    }

    output{
        Array[File] tar_gz_pe = runIGV_whole_genome.pe_plots
    }
}

task runIGV_whole_genome{
    input{
        Array[File] varfiles
        File fasta
        File fasta_idx
        File fasta_dict
        File nested_repeats
        File simple_repeats
        File empty_track
        File sample_crai_cram
        Array[String] families
        File ped_file
        Array[Array[String]] samples
        Array[Array[File]] crams
        Array[Array[File]] crais
        String buffer
        String buffer_large
        String igv_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([varfiles, fasta, fasta_idx, fasta_dict, nested_repeats, simple_repeats, empty_track, ped_file, crams, crais]), "GB")
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
            #export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
            mkdir pe_igv_plots
            i=0
            for varfile in ~{sep=' ' varfiles}
            do
                let "i=$i+1"
                family = ${families[$i]}
                samples_file = write_lines(${samples[$i]})
                crams_file = write_lines(${crams[$i]})
                crai_list = ${crais[$i]}
                #if HB does not think above will work we can do
                #grep -w family ~{ped_file} | cut -f1 | sort -u > samples.txt
                #grep -w -f samples.txt ~{sample_crai_cram} | cut -f3 > crams.txt
                python /src/makeigvpesr.py -v "${varfile}" -n ~{nested_repeats} -s ~{simple_repeats} -e ~{empty_track} -f ~{fasta} -fam_id ${family} -samples $samples_file -crams $crams_file -p ~{ped_file} -o pe_igv_plots -b ~{buffer} -l ~{buffer_large}
                bash pe.sh
                xvfb-run --server-args="-screen 0, 1920x540x24" bash /IGV_2.4.14/igv.sh -b pe.$i.txt
                tar -czf $family_pe_igv_plots.tar.gz pe_igv_plots
            done;

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
        Array[File] pe_plots="$family_pe_igv_plots.tar.gz"
        Array[File] pe_txt = "pe.$i.txt"
        }
    }