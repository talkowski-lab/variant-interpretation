#Author: Alba Sanchis-Juan <asanchis@broadinstitute.org>
#This script makes SV plots from PE SR and RD evidence


version 1.0

import "Structs.wdl"

workflow IGV_evidence {
    input{
        File varfile
        String family
        Array[String] samples
        Array[File] disc_files
        Array[File] split_files
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_reformat_pe
        RuntimeAttr? runtime_attr_reformat_sr
        RuntimeAttr? runtime_attr_update_pe_sr
    }

    scatter (i in range(length(samples))){
        call reformatPE{
            input:
                varfile=varfile,
                sample = samples[i],
                disc_file=disc_files[i],
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override=runtime_attr_reformat_pe

        }

        call reformatSR{
            input:
                varfile=varfile,
                sample = samples[i],
                split_file=split_files[i],
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override=runtime_attr_reformat_sr
        }
    }

    call update_sample_pe_sr{
        input:
            pe = reformatPE.pe_reformat_string,
            sr = reformatSR.sr_reformat_string,
            samples = samples,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_override = runtime_attr_update_pe_sr
        }

    output{
        Array[File] pe_files = reformatPE.pe_reformat
        Array[File] sr_files = reformatSR.sr_reformat
        Array[String] pe_strings = reformatPE.pe_reformat_string
        Array[String] sr_strings = reformatSR.sr_reformat_string
        File updated_sample_pe_sr = update_sample_pe_sr.changed_sample_pe_sr   
    }
}

task reformatPE{
    input{
        File varfile
        String sample
        File disc_file
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([varfile, disc_file]), "GB")
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

            cat ~{varfile} | cut -f1-3 | awk '{if ($3-$2>=15000) print $1"\t"$2"\t"$2 "\n" $1"\t"$3"\t"$3;else print}' | awk '{$2-=3000}1' OFS='\t' | awk '{$3+=3000}1' OFS='\t' | sort -k1,1 -k2,2n | bgzip -c > regions.bed.gz
            tabix -p bed regions.bed.gz

            ##Reformat disc reads
            zcat ~{disc_file} | \
                awk '$1 == $4 {print $1"\t"$2"\t"$5"\t.\t1\t"$3"\t255,0,0"}' | \
                sed -e 's/-\t255,0,0/-\t0,255,0/g' | \
                bgzip -c >  tmp.bed.gz

            tabix -p bed tmp.bed.gz
            bedtools intersect -a tmp.bed.gz -b regions.bed.gz -f 0.5 -wa > ~{sample}.pe_roi.junctions.bed
            >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
    output{
        File pe_reformat="~{sample}.pe_roi.junctions.bed"
        String pe_reformat_string="~{sample}.pe_roi.junctions.bed"
    }
}


task reformatSR{
    input{
        File varfile
        String sample
        File split_file
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([varfile, split_file]), "GB")
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

            cat ~{varfile} | cut -f1-3 | awk '{if ($3-$2>=15000) print $1"\t"$2"\t"$2 "\n" $1"\t"$3"\t"$3;else print}' | awk '{$2-=3000}1' OFS='\t' | awk '{$3+=3000}1' OFS='\t' | sort -k1,1 -k2,2n | bgzip -c > regions.bed.gz
            tabix -p bed regions.bed.gz

            ##Reformat split reads
            zcat ~{split_file} | awk -F"\t" 'BEGIN { OFS="\t" }{print $1,$2-5,$2+5,$3,$4}' | \
                sed -e 's/left/-/g' | \
                sed -e 's/right/+/g' | \
                bgzip -c > tmp.bed.gz

            tabix -p bed tmp.bed.gz
            tabix -R regions.bed.gz tmp.bed.gz | \
                bgzip -c > ~{sample}.sr_roi.junctions.bed
            >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
    output{
        File sr_reformat="~{sample}.sr_roi.junctions.bed"
        String sr_reformat_string="~{sample}.sr_roi.junctions.bed"
    }
}

task update_sample_pe_sr{
    input {
        Array[String] pe
        Array[String] sr
        Array[String] samples
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }
    
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
           pe = ~{write_lines(pe)}
           sr = ~{write_lines(sr)}
           samples = ~{write_lines(samples)}
           paste samples pe sr > updated_samples_pe_sr.txt

        >>>

    output{
        File changed_sample_pe_sr = "updated_samples_pe_sr.txt"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: variant_interpretation_docker
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}
