version 1.0
    
import "https://raw.githubusercontent.com/broadinstitute/gatk-sv/v0.26-beta/wdl/Structs.wdl"

workflow fixCountsCoords {

    input {
        File counts_fof
        String docker
        RuntimeAttr? runtime_attr_override

    }

    Array[String] counts_files = transpose(read_tsv(counts_fof))[0]

    scatter (file in counts_files) {
        call fixCounts {
            input:
                counts_file=file,
                docker=docker
        }
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
        preemptible: 3,
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
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}
