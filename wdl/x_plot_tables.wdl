version 1.0

import "https://raw.githubusercontent.com/broadinstitute/gatk-sv/v0.26-beta/wdl/Structs.wdl"


workflow x_plot_tables {

    input {

        File ped_input
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
        Array[String] sample_list
        Array[String] cram_files
    }   

    call x_makeDataTable {
        input:
            ped_input = ped_input,
            sample_list = sample_list,
            cram_files = cram_files,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_override
    }
}


task x_makeDataTable{
    input{
        Array[String] sample_list
        Array[String] cram_files
        File ped_input
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
        File output_table = "family.sample.table.tsv"
    }

    File samples_list = write_lines(sample_list)
    File cram_list = write_lines(cram_files)

    command {

        Rscript /src/variant-interpretation/scripts/create_familyid_table.R ${samples_list} ${cram_list} ${ped_input} family.sample.table.tsv

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

