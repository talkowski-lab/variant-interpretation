version 1.0
    
import "https://raw.githubusercontent.com/broadinstitute/gatk-sv/v0.26-beta/wdl/Structs.wdl"

workflow fixCountsCoords {

    input {
        File in_file
        String docker
        RuntimeAttr? runtime_attr_override
    }

    call fixCounts {
        input:
            counts_file=in_file,
            docker=docker
    }

    output {
        File output_file_name=fixCounts.output_file
    }

}

task fixCounts{
    input{
        File counts_file
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 12,
        boot_disk_gb: 8,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String output_name=sub(basename(counts_file), "txt.gz", "coords.txt.gz")

    output{
        File output_file=output_name
    }

    command <<<

        zcat ~{counts_file} | \
            awk -v OFS='\t' -F'\t' '/^[^@|CONTIG]/{ $2 = $2+1 }{print $0}' | \
            bgzip -c > ~{output_name}
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}
