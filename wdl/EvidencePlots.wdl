#Author: Alba Sanchis-Juan <asanchis@broadinstitute.org>
#This script makes SV plots from PE SR and RD evidence


version 1.0

import "Structs.wdl"

workflow IGV_evidence {
    input{
        File varfile
        Array[String] samples
        Array[File] disc_files
        Array[File] split_files
        Boolean file_localization
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_reformat_pe
        RuntimeAttr? runtime_attr_reformat_sr
        RuntimeAttr? runtime_attr_update_pe_sr
    }

    scatter (i in range(length(samples))){
        if(!(file_localization)){
            call reformatPE_parse{
                input:
                    varfile=varfile,
                    sample = samples[i],
                    disc_file=disc_files[i],
                    variant_interpretation_docker=variant_interpretation_docker,
                    runtime_attr_override=runtime_attr_reformat_pe

            }

            call reformatSR_parse{
                input:
                    varfile=varfile,
                    sample = samples[i],
                    split_file=split_files[i],
                    variant_interpretation_docker=variant_interpretation_docker,
                    runtime_attr_override=runtime_attr_reformat_sr
            }
        }

        if (file_localization){
            call reformatPE_localize{
                input:
                    varfile=varfile,
                    sample = samples[i],
                    disc_file=disc_files[i],
                    variant_interpretation_docker=variant_interpretation_docker,
                    runtime_attr_override=runtime_attr_reformat_pe

            }

            call reformatSR_localize{
                input:
                    varfile=varfile,
                    sample = samples[i],
                    split_file=split_files[i],
                    variant_interpretation_docker=variant_interpretation_docker,
                    runtime_attr_override=runtime_attr_reformat_sr
            }
        }
    }

    call update_sample_pe_sr{
        input:
            pe = select_first([reformatPE_parse.pe_reformat, reformatPE_localize.pe_reformat]),
            sr = select_first([reformatSR_parse.sr_reformat, reformatSR_localize.sr_reformat]),
            samples = samples,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_override = runtime_attr_update_pe_sr
        }

    output{
        Array[File] pe_files = select_first([reformatPE_parse.pe_reformat, reformatPE_localize.pe_reformat])
        Array[File] sr_files = select_first([reformatSR_parse.sr_reformat, reformatSR_localize.sr_reformat])
        Array[String] pe_strings = select_first([reformatPE_parse.pe_reformat_string, reformatPE_localize.pe_reformat_string])
        Array[String] sr_strings = select_first([reformatSR_parse.sr_reformat_string, reformatSR_localize.sr_reformat_string])
        File updated_sample_pe_sr = update_sample_pe_sr.changed_sample_pe_sr   
    }
}

task reformatPE_parse{
    input{
        File varfile
        String sample
        String disc_file
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([varfile]), "GB")
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

            export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
            tabix -R regions.bed.gz ~{disc_file} | bgzip -c > subset_disc.txt.gz  
            tabix -p bed subset_disc.txt.gz
            ##Reformat disc reads
            zcat subset_disc.txt.gz | \
                grep -w ~{sample} | \
                awk '$1 == $4 {print $1"\t"$2"\t"$5"\t.\t1\t"$3"\t255,0,0"}' | \
                sed -e 's/-\t255,0,0/-\t0,255,0/g' > ~{sample}.pe_roi.junctions.bed
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


task reformatSR_parse{
    input{
        File varfile
        String sample
        String split_file
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([varfile]), "GB")
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

            export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
            tabix -R regions.bed.gz ~{split_file} | bgzip -c > subset_sr.txt.gz  
            tabix -p bed subset_sr.txt.gz
            ##Reformat split reads
            zcat subset_sr.txt.gz | grep -w ~{sample} | \
            awk -F"\t" 'BEGIN { OFS="\t" }{print $1,$2-5,$2+5,$3,$4}' | \
                sed -e 's/left/-/g' | \
                sed -e 's/right/+/g' | \
                bgzip -c > ~{sample}.sr_roi.bed.gz

            tabix -p bed ~{sample}.sr_roi.bed.gz
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
        File sr_reformat="~{sample}.sr_roi.bed.gz"
        String sr_reformat_string="~{sample}.sr_roi.bed.gz"
    }
}

task reformatPE_localize{
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
                grep -w ~{sample} | \
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


task reformatSR_localize{
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
            zcat ~{split_file} | grep -w ~{sample} | \
            awk -F"\t" 'BEGIN { OFS="\t" }{print $1,$2-5,$2+5,$3,$4}' | \
                sed -e 's/left/-/g' | \
                sed -e 's/right/+/g' | \
                bgzip -c > tmp.bed.gz

            tabix -p bed tmp.bed.gz
            tabix -R regions.bed.gz tmp.bed.gz | \
                bgzip -c > ~{sample}.sr_roi.bed.gz
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
        File sr_reformat="~{sample}.sr_roi.bed.gz"
        String sr_reformat_string="~{sample}.sr_roi.bed.gz"
    }
}

task update_sample_pe_sr{
    input {
        Array[File] pe
        Array[File] sr
        Array[String] samples
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
      pe: {
        localization_optional: true
      }
      sr: {
        localization_optional: true
      }

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
        
           paste ~{write_lines(samples)} ~{write_lines(pe)} ~{write_lines(sr)} > updated_samples_pe_sr.txt

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
