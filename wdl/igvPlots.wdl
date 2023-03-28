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
        Array[String] crams
        Array[String] crais
        String buffer
        String buffer_large
        String igv_docker
        RuntimeAttr? runtime_attr_igv
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
            buffer_large = buffer_large,
            igv_docker = igv_docker,
            runtime_attr_override = runtime_attr_igv
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
        Array[String] crams
        Array[String] crais
        String buffer
        String buffer_large
        String igv_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([varfile, fasta, fasta_idx, fasta_dict, nested_repeats, simple_repeats, empty_track, ped_file, crams, crais]), "GB")
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
            mkdir pe_igv_plots
            cat ~{varfile} | cut -f1-3 | awk '{$2-=3000}1' OFS='\t' | awk '{$3+=3000}1' OFS='\t' > regions.bed
            #localize cram files
            for cram in ~{sep=' ' crams}
            do
                export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
                samtools view -h -o new.$cram $cram -L regions.bed
                samtools index new.$cram
            done
            ls new.* > crams.txt
        python /src/makeigvpesr.py -v ~{varfile} -n ~{nested_repeats} -s ~{simple_repeats} -e ~{empty_track} -f ~{fasta} -fam_id ~{family} -samples ~{sep="," samples} -crams crams.txt -p ~{ped_file} -o pe_igv_plots -b ~{buffer} -l ~{buffer_large} -i pe.txt -bam pe.sh
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
        File pe_plots="~{family}_pe_igv_plots.tar.gz"
        File pe_txt = "pe.txt"
        }
    }