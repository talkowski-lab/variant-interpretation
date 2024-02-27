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
        Array[File]? vep_vcf_files
        Array[File]? vep_annotated_final_vcf
        File reference_tsv
        File full_input_tsv
        File hg38_reference
        File hg38_reference_dict
        File hg38_reference_fai
        Float snv_scale=1
        Float indel_scale=1
        String jvarkit_docker
        String vep_hail_docker
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
            scale=snv_scale,
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
            scale=indel_scale,
            num_variants=getNumIndels.num_variants
    }

    # only annotate Indels for now
    Array[File] vep_files = select_first([vep_vcf_files, vep_annotated_final_vcf])
    call convertTSVtoVCF {
        input:
        tsv=downsampleIndels.downsampled_tsv,
        vcf_file=vep_files[0],
        vep_hail_docker=vep_hail_docker
    }

    call annotatePOLYX {
        input:
        vcf_file=convertTSVtoVCF.output_vcf,
        hg38_reference=hg38_reference,
        hg38_reference_fai=hg38_reference_fai,
        hg38_reference_dict=hg38_reference_dict,
        jvarkit_docker=jvarkit_docker
    }

    output {
        File downsampled_tsv_SNV = downsampleSNVs.downsampled_tsv
        File downsampled_tsv_Indel = downsampleIndels.downsampled_tsv
        File downsampled_polyx_vcf_Indel = annotatePOLYX.polyx_vcf
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
        Float scale
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

    Int desired_num_variants = ceil(num_variants * scale)

    command <<<
        type_col=$(cat ~{full_input_tsv} | head -n 1 | awk -v name='TYPE' '{for (i=1;i<=NF;i++) if ($i==name) print i; exit}')
        var_type_tsv=$(basename ~{full_input_tsv} '.tsv')_~{var_type}.tsv
        cat ~{full_input_tsv} | awk -v col=$type_col '{ if ($col == "~{var_type}") { print } }' > $var_type_tsv
        tot_n_vars=$(cat $var_type_tsv | wc -l)
        echo "total $tot_n_vars ~{var_type}s"
        desired_num_variants=~{desired_num_variants}
        num_downsample=$(( desired_num_variants < tot_n_vars ? desired_num_variants : tot_n_vars ))
        percent_downsample=$(echo "scale=4 ; $num_downsample / $tot_n_vars" | bc)
        echo "downsampling $percent_downsample"
        downsampled_tsv=$(basename $var_type_tsv '.tsv')_downsampled.tsv
        cat ~{full_input_tsv} | head -n 1 > $downsampled_tsv
        awk -v percent_downsample=$percent_downsample 'BEGIN {srand()} !/^$/ { if (rand() <= percent_downsample) print $0}' $var_type_tsv >> $downsampled_tsv
    >>>

    output {
        File downsampled_tsv = "~{basename(full_input_tsv, '.tsv')}_~{var_type}_downsampled.tsv"
    }
}


task convertTSVtoVCF {
    input {
        File tsv
        File vcf_file  # for header
        String vep_hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(tsv, "GB") + size(vcf_file, "GB")
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
    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])

    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: vep_hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        cat <<EOF > tsv_to_vcf.py 
        import os
        import sys
        import pandas as pd
        import numpy as np
        import hail as hl

        tsv = sys.argv[1]
        vcf_file = sys.argv[2]
        cores = sys.argv[3]
        mem = int(np.floor(float(sys.argv[4])))

        hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

        df = pd.read_csv(tsv, sep='\t')[['CHROM','POS','REF','ALT']].to_csv('temp.tsv', sep='\t',index=False)
        mt = hl.import_matrix_table('temp.tsv', row_fields={'CHROM':'str','POS':'int','REF':'str','ALT':'str'})

        mt = mt.annotate_rows(locus=hl.locus(mt.CHROM, mt.POS, 'GRCh38'),
                        alleles=hl.array([mt.REF, mt.ALT]))
        mt = mt.key_rows_by('locus','alleles')

        header = hl.get_vcf_metadata(vcf_file)
        hl.export_vcf(mt, os.path.basename(tsv).split('.tsv')[0]+'.vcf', metadata=header)
        EOF

        python3.9 tsv_to_vcf.py ~{tsv} ~{vcf_file} ~{cpu_cores} ~{memory}
    >>>

    output {
        File output_vcf = basename(tsv, '.tsv') + '.vcf'
    }
}

task annotatePOLYX {
    input {
        File vcf_file
        File hg38_reference
        File hg38_reference_fai
        File hg38_reference_dict
        String jvarkit_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 2.0
    Float input_mem_scale = 3.0
    Float input_disk_scale = 5.0
    
    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb + input_size * input_mem_scale,
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
        docker: jvarkit_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String out_vcf = basename(vcf_file, '.vcf')+'_POLYX.vcf'

    command {
        java -jar /opt/jvarkit/dist/jvarkit.jar vcfpolyx -R ~{hg38_reference} -o ~{out_vcf} ~{vcf_file}
    }

    output {
        File polyx_vcf = out_vcf
    }
}