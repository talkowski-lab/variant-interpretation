version 1.0

import "Structs.wdl"

workflow createBaiOrCrai {
    input {

        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
        File bam_or_cram_file
        
    }   

    if (basename(bam_or_cram_file, ".cram") + ".cram" == basename(bam_or_cram_file)) {
        call x_createCrais{
            input:
                cram_file = bam_or_cram_file,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_override
        }
    }
   
    if (basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)) {
        call x_createBais{
            input:
                bam_file = bam_or_cram_file,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_override
        }
    }

    output{
        File bam_or_cram_index = select_first([x_createCrais.crai_file, x_createBais.bai_file])
    }
}

    task x_createCrais{
    input{
        File cram_file
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }
    RuntimeAttr default_attr = object {
        cpu: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible: 3,
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
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task x_createBais{
    input{
        File bam_file
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }
    RuntimeAttr default_attr = object {
        cpu: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible: 3,
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
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}
