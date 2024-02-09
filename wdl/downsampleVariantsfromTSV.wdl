version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow downsampleVariantsfromTSV {
    input {
        File reference_tsv
        File full_input_tsv
        String sv_base_mini_docker
    }

    call getNumVars as getNumSNVs {
        input:
            reference_tsv=reference_tsv,
            sv_base_mini_docker=sv_base_mini_docker,
            var_type='SNV'
    }
    call downsampleVariants as downsampleSNVs {
        input:
            full_input_tsv=full_input_tsv,
            sv_base_mini_docker=sv_base_mini_docker,
            var_type='SNV',
            num_variants=getNumSNVs.num_variants
    }

    call getNumVars as getNumIndels {
        input:
            reference_tsv=reference_tsv,
            sv_base_mini_docker=sv_base_mini_docker,
            var_type='Indel'
    }
    call downsampleVariants as downsampleIndels {
        input:
            full_input_tsv=full_input_tsv,
            sv_base_mini_docker=sv_base_mini_docker,
            var_type='Indel',
            num_variants=getNumIndels.num_variants
    }

    output {
        File downsampled_tsv_SNV = downsampleSNVs.downsampled_tsv
        File downsampled_tsv_Indel = downsampleIndels.downsampled_tsv
    }
}

task getNumVars {
    input {
        File reference_tsv
        String sv_base_mini_docker
        String var_type
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(reference_tsv, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        type_col=$(cat ~{reference_tsv} | head -n 1 | awk -v name='TYPE' '{for (i=1;i<=NF;i++) if ($i==name) print i; exit}')
        n_vars=$(cat ~{reference_tsv} | awk -v col=$type_col '{ if ($col == "~{var_type}") { print } }' | wc -l)
        echo $n_vars > n_vars.txt
    >>>

    output {
        Int num_variants = read_lines('n_vars.txt')[0]
    }
}

task downsampleVariants {
    input {
        File full_input_tsv
        String sv_base_mini_docker
        String var_type
        Int num_variants
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(full_input_tsv, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -eou pipefail
        type_col=$(cat ~{full_input_tsv} | head -n 1 | awk -v name='TYPE' '{for (i=1;i<=NF;i++) if ($i==name) print i; exit}')
        var_type_tsv=$(basename ~{full_input_tsv} '.tsv')_~{var_type}.tsv
        cat ~{full_input_tsv} | awk -v col=$type_col '{ if ($col == "~{var_type}") { print } }' > $var_type_tsv
        tot_n_vars=$(cat $var_type_tsv | wc -l)
        echo "total $tot_n_vars ~{var_type}s"
        percent_downsample=$(echo "scale=4 ; ~{num_variants} / $tot_n_vars" | bc)
        echo "downsampling $percent_downsample"
        downsampled_tsv=$(basename $var_type_tsv '.tsv')_downsampled.tsv
        cat ~{full_input_tsv} | head -n 1 > $downsampled_tsv
        awk -v percent_downsample=$percent_downsample 'BEGIN {srand()} !/^$/ { if (rand() <= percent_downsample) print $0}' $var_type_tsv >> $downsampled_tsv
    >>>

    output {
        File downsampled_tsv = "~{basename(full_input_tsv, '.tsv')}_~{var_type}_downsampled.tsv"
    }
}