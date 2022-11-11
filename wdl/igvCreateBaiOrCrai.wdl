version 1.0

import "https://raw.githubusercontent.com/broadinstitute/gatk-sv/v0.26-beta/wdl/Structs.wdl"

workflow createBaiOrCrai {
    input {

        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
        File bam_or_cram_file
        
    }   

    if (basename(bam_or_cram_file, ".cram") + ".cram" == basename(bam_or_cram_file)) {
        call createCrais{
            input:
                cram_file = bam_or_cram_file,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_override
        }
    }
   
    if (basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)) {
        call createBais{
            input:
                bam_file = bam_or_cram_file,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_override
        }
    }

    output{
        File bam_or_cram_index = select_first([createCrais.crai_file, createBais.bai_file])
    }
}

    task createCrais{
    input{
        File cram_file
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible_tries: 3,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    output{
        File crai_file = "${cram_file}.crai"
    }
    command {
        samtools index -b ${cram_file}
    }
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task createBais{
    input{
        File bam_file
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible_tries: 3,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    output{
        File bai_file = "${bam_file}.bai"
    }
    command {
        samtools index -b ${bam_file}
    }
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}
