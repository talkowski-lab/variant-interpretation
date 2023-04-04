version 1.0

##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

import "Structs.wdl"

workflow IGV {
    input{
        File varfile
        String family
        File ped_file
        Array[String] samples
        Array[String] crams
        Array[String] crais
        File reference
        String buffer
        String buffer_large
        String igv_docker
        RuntimeAttr? runtime_attr_igv
    }

    call runIGV_whole_genome{
        input:
            varfile = varfile,
            family = family,
            ped_file = ped_file,
            samples = samples,
            crams = crams,
            crais = crais,
            buffer = buffer,
            buffer_large = buffer_large,
            reference = reference,
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
        File reference
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

    Float input_size = size(select_all([varfile, reference, ped_file]), "GB")
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
            cat ~{varfile} | cut -f1-3 | awk '{if ($3-$2>=15000) print $1"\t"$2"\t"$2 "\n" $1"\t"$3"\t"$3;else print}' | awk '{$2-=3000}1' OFS='\t' | awk '{$3+=3000}1' OFS='\t' | sort -k1,1 -k2,2n | bgzip -c > regions.bed.gz
            tabix -p bed regions.bed.gz
            #localize cram files
            for cram in ~{sep=' ' crams}
            do
                name=$(echo $cram|awk -F"/" '{print $NF}'|sed 's/.cram//g')
                export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
                samtools view -h -C -T ~{reference} -o new.$name.cram $cram -L regions.bed.gz -M
                samtools index new.$name.cram
            done
            ls *.cram > crams.txt

            i=0
            while read -r line
            do
                let "i=$i+1"
                echo "$line" > new.varfile.$i.bed
                python /src/makeigvpesr.py -v new.varfile.$i.bed -fam_id ~{family} -samples ~{sep="," samples} -crams crams.txt -p ~{ped_file} -o pe_igv_plots -b ~{buffer} -l ~{buffer_large} -i pe.$i.txt -bam pe.$i.sh
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
        Array[File] crams = glob("*.cram")
        Array[File] crais = glob("*.crai")
        Array[File] varfile = glob("new.varfile.*.bed")
        }
    }