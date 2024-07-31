version 1.0

##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

import "Structs2.wdl"

workflow IGV {
    input{
        File varfile
        String family
        File ped_file
        Array[String] samples
        Boolean file_localization
        Boolean requester_pays
        Array[File]? crams_localize
        Array[File]? crais_localize
        File reference
        File reference_index
        Int igv_max_window
        String buffer
        String igv_docker
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_igv
        RuntimeAttr? runtime_attr_localize_reads
        Array[String]? crams_parse
        Array[String]? crais_parse
        File? updated_sample_crai_cram
        File? sample_crai_cram
        
    }

    if (file_localization) {
        Array[File] crams_localize_ = select_first([crams_localize])
        Array[File] crais_localize_ = select_first([crais_localize])
        File sample_crai_cram_ = select_first([sample_crai_cram])

        if (requester_pays){
        # move the reads nearby -- handles requester_pays and makes cross-region transfers just once
            scatter(i in range(length(crams_localize_))) {
                call LocalizeReads as LocalizeReadsLocalize{
                    input:
                        reads_path = crams_localize_[i],
                        reads_index = crais_localize_[i],
                        runtime_attr_override = runtime_attr_localize_reads
                }
            }
        }

        call runIGV_whole_genome_localize {
            input:
                varfile = varfile,
                family = family,
                ped_file = ped_file,
                samples = samples,
                crams = select_first([LocalizeReadsLocalize.output_file, crams_localize_]),
                crais = select_first([LocalizeReadsLocalize.output_index, crais_localize_]),
                sample_crai_cram = sample_crai_cram_,
                buffer = buffer,
                reference = reference,
                reference_index = reference_index,
                igv_max_window = igv_max_window,
                igv_docker = igv_docker,
                runtime_attr_override = runtime_attr_igv
        }
    }

    if (!(file_localization)) {
        Array[String] crams_parse_ = select_first([crams_parse])
        Array[String] crais_parse_ = select_first([crais_parse])
        File updated_sample_crai_cram_ = select_first([updated_sample_crai_cram])

        call runIGV_whole_genome_parse{
            input:
                varfile = varfile,
                family = family,
                ped_file = ped_file,
                samples = samples,
                updated_sample_crai_cram = updated_sample_crai_cram_,
                buffer = buffer,
                reference = reference,
                reference_index = reference_index,
                igv_max_window = igv_max_window,
                igv_docker = igv_docker,
                runtime_attr_override = runtime_attr_igv
        }
    }

    output{
        File tar_gz_pe = select_first([runIGV_whole_genome_localize.pe_plots, runIGV_whole_genome_parse.pe_plots])
    }
}

task runIGV_whole_genome_localize{
        input{
            File varfile
            File reference
            File reference_index
            Int igv_max_window
            String family
            File ped_file
            Array[String] samples
            Array[File] crams
            Array[File] crais
            File sample_crai_cram
            String buffer
            String igv_docker
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
            mkdir pe_igv_plots
            head -n+1 ~{ped_file} > family_ped.txt
            grep -w ~{family} ~{ped_file} >> family_ped.txt
            python3.6 /src/renameCramsLocalize.py --ped family_ped.txt --scc ~{sample_crai_cram}
            cut -f4 changed_sample_crai_cram.txt > crams.txt
            
            while read sample crai cram new_cram new_crai
            do
                mv $cram $new_cram
                mv $crai $new_crai
            done<changed_sample_crai_cram.txt

            i=0
            while read -r line
            do
                let "i=$i+1"
                echo "$line" > new.varfile.$i.bed
                python /src/variant-interpretation/scripts/makeigvpesr.py -v new.varfile.$i.bed -fam_id ~{family} -samples ~{sep="," samples} -crams crams.txt -p ~{ped_file} -o pe_igv_plots -b ~{buffer}  -i pe.$i.txt -bam pe.$i.sh -m ~{igv_max_window}
                bash pe.$i.sh
                xvfb-run --server-args="-screen 0, 1920x540x24" bash /IGV_Linux_2.16.0/igv.sh -b pe.$i.txt
            done < ~{varfile}
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
        Array[File] pe_txt = glob("pe.*.txt")
        Array[File] pe_sh = glob("pe.*.sh")
        Array[File] varfile = glob("new.varfile.*.bed")
        }
    }

task runIGV_whole_genome_parse{
    input{
        File varfile
        File reference
        File reference_index
        Int igv_max_window
        String family
        File ped_file
        Array[String] samples
        File updated_sample_crai_cram
        String buffer
        String igv_docker
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
            mkdir pe_igv_plots
            cat ~{varfile} | cut -f1-3 | awk '{if (($3-$2)+int(($3-$2)*1.5)>=~{igv_max_window}) print $1"\t"$2-~{buffer}"\t"$2+~{buffer} "\n" $1"\t"$3-~{buffer}"\t"$3+~{buffer};
                else print $1"\t"($2-int(($3-$2)*0.25))-~{buffer}"\t"$3+int(($3-$2)*0.25)+~{buffer}}' | sort -k1,1 -k2,2n | bgzip -c > regions.bed.gz
            tabix -p bed regions.bed.gz
            #localize cram files
            while read sample crai cram new_cram new_crai
            do
                name=$( echo $new_cram|awk -F"/" '{print $NF}'|sed 's/.cram//g' )
                gsutil cp $crai $( basename $cram | sed 's/\.cram$/.crai/g' )
                export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
                samtools view -h -C -T ~{reference} -o $name.cram $cram -L regions.bed.gz -M
                samtools index $name.cram
            done<~{updated_sample_crai_cram}
            ls *.cram > crams.txt

            i=0
            while read -r line
            do
                let "i=$i+1"
                echo "$line" > new.varfile.$i.bed
                python /src/variant-interpretation/scripts/makeigvpesr.py -v new.varfile.$i.bed -fam_id ~{family} -samples ~{sep="," samples} -crams crams.txt -p ~{ped_file} -o pe_igv_plots -b ~{buffer} -i pe.$i.txt -bam pe.$i.sh -m ~{igv_max_window}
                bash pe.$i.sh
                xvfb-run --server-args="-screen 0, 1920x540x24" bash /IGV_Linux_2.16.0/igv.sh -b pe.$i.txt
            done < ~{varfile}
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
        Array[File] pe_txt = glob("pe.*.txt")
        Array[File] pe_sh = glob("pe.*.sh")
        Array[File] varfile = glob("new.varfile.*.bed")
        Array[File] cramfiles = glob("*cram")
        }
    }

task LocalizeReads {
  input {
    File reads_path
    File reads_index
    RuntimeAttr? runtime_attr_override
  }

#  Float input_size = size(reads_path, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
#                                  disk_gb: ceil(50.0 + input_size),
                                  disk_gb: ceil(60.0),
                                  cpu: 2,
                                  preemptible: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
                                
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu, runtime_default.cpu])
    preemptible: select_first([runtime_override.preemptible, runtime_default.preemptible])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: "ubuntu:18.04"
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

#  Int disk_size = ceil(50 + size(reads_path, "GB"))

  command {
      set -exuo pipefail

      cp ~{reads_path} $(basename ~{reads_path})
      cp ~{reads_index} $(basename ~{reads_index})
#    ln -s ~{reads_path}
#    ln -s ~{reads_index}
  }
  output {
    File output_file = basename(reads_path)
    File output_index = basename(reads_index)
  }
}