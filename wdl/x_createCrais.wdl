version 1.0

import "Structs.wdl"

workflow createCrais {
    input {

        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
        File cram_file
        
    }   
    
    call x_createCrais{
        input:
            cram_file = cram_file,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_override
    }

    output{
        File cram_index = x_createCrais.crai_file
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
