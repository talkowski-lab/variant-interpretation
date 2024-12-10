version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow removeRegionsVariantsfromTSV {
    input {
        File full_input_tsv
        File remove_regions_bed
        File hg38_reference
        File hg38_reference_dict
        File hg38_reference_fai
        String jvarkit_docker
        String hail_docker
        String genome_build='GRCh38'
        Boolean prioritize_coding=true
        RuntimeAttr? runtime_attr_remove_regions
    }


    scatter (var_type in ['SNV', 'Indel']) {

        call removeRegionsVariants {
            input:
            remove_regions_bed=remove_regions_bed,
            full_input_tsv=full_input_tsv,
            hail_docker=hail_docker,
            var_type=var_type,
            genome_build=genome_build,
            runtime_attr_override=runtime_attr_remove_regions,
            prioritize_coding=prioritize_coding
        }

        call convertTSVtoVCF {
            input:
            tsv=removeRegionsVariants.filtered_tsv,
            runtime_attr_override=runtime_attr_remove_regions,
            hail_docker=hail_docker,
            genome_build=genome_build
        }

        call annotatePOLYX {
            input:
            vcf_file=convertTSVtoVCF.output_vcf,
            hg38_reference=hg38_reference,
            hg38_reference_fai=hg38_reference_fai,
            hg38_reference_dict=hg38_reference_dict,
            runtime_attr_override=runtime_attr_remove_regions,
            jvarkit_docker=jvarkit_docker
        }

        call mergePOLYX {
            input:
            polyx_vcf=annotatePOLYX.polyx_vcf,
            tsv=removeRegionsVariants.filtered_tsv,
            runtime_attr_override=runtime_attr_remove_regions,
            hail_docker=hail_docker
        }
    }

    output {
        File filtered_tsv_SNV = mergePOLYX.filtered_polyx_tsv[0]
        File filtered_tsv_Indel = mergePOLYX.filtered_polyx_tsv[1]
    }
}

task convertTSVtoVCF {
    input {
        File tsv
        String hail_docker
        String genome_build
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(tsv, "GB") 
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
        docker: hail_docker
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
        cores = sys.argv[2]
        mem = int(np.floor(float(sys.argv[3])))
        build = sys.argv[4]

        hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

        df = pd.read_csv(tsv, sep='\t')[['CHROM','POS','REF','ALT']].to_csv('temp.tsv', sep='\t',index=False)
        mt = hl.import_matrix_table('temp.tsv', row_fields={'CHROM':'str','POS':'int','REF':'str','ALT':'str'})

        mt = mt.annotate_rows(locus=hl.locus(mt.CHROM, mt.POS, build),
                        alleles=hl.array([mt.REF, mt.ALT]))
        mt = mt.key_rows_by('locus','alleles')

        hl.export_vcf(mt, os.path.basename(tsv).split('.tsv')[0]+'.vcf')
        EOF

        python3 tsv_to_vcf.py ~{tsv} ~{cpu_cores} ~{memory} ~{genome_build}
    >>>

    String file_ext = if sub(basename(tsv), '\\.gz', '')==basename(tsv) then '.tsv' else '.tsv.gz'
    output {
        File output_vcf = basename(tsv, file_ext) + '.vcf'
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

task removeRegionsVariants {
    input {
        File full_input_tsv
        File remove_regions_bed
        String var_type
        String hail_docker
        String genome_build
        Boolean prioritize_coding
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
    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String file_ext = if sub(basename(full_input_tsv), '\\.gz', '')==basename(full_input_tsv) then '.tsv' else '.tsv.gz'
    String output_name = "~{basename(full_input_tsv, file_ext)}_~{var_type}.removed.regions.tsv.gz"

    command <<<
        cat <<EOF > removeRegions.py
        import pandas as pd
        import numpy as np
        import os
        import sys
        import ast
        import hail as hl

        full_input_tsv = sys.argv[1]
        var_type = sys.argv[2]
        remove_regions_bed = sys.argv[3]
        output_name = sys.argv[4]
        prioritize_coding = ast.literal_eval(sys.argv[5].capitalize())
        cores = sys.argv[6]
        mem = int(np.floor(float(sys.argv[7])))
        build = sys.argv[8]

        hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

        ht = hl.import_table(full_input_tsv, force=full_input_tsv.split('.')[-1] in ['gz', 'bgz'])
        exclude_ht = hl.import_bed(remove_regions_bed, reference_genome=build)

        ht = ht.annotate(locus=hl.parse_variant(ht.ID, build).locus,
                        alleles=hl.parse_variant(ht.ID, build).alleles)

        ht = ht.filter(ht.TYPE==var_type)

        if prioritize_coding:  # keep all coding
            ht = ht.filter((~hl.is_defined(exclude_ht[ht.locus])) | (hl.bool(ht.isCoding)))
        else:
            ht = ht.filter(hl.is_defined(exclude_ht[ht.locus]), keep=False)
        ht.export(output_name)
        EOF

        python3 removeRegions.py ~{full_input_tsv} ~{var_type} ~{remove_regions_bed} ~{output_name} ~{prioritize_coding} \
            ~{cpu_cores} ~{memory} ~{genome_build}
    >>>

    output {
        File filtered_tsv = output_name
    }
}

task mergePOLYX {
    input {
        File polyx_vcf
        File tsv
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size([polyx_vcf, tsv], "GB")
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
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        cat <<EOF > mergePOLYX.py
        import pandas as pd
        import os
        import sys

        polyx_vcf = sys.argv[1]
        tsv = sys.argv[2]

        df = pd.read_csv(tsv, sep='\t')
        polyx_df = pd.read_csv(polyx_vcf, sep='\t', comment='#', header=None)
        df['POLYX'] = polyx_df[7].str.split('=').str[-1].astype(int)
        df.to_csv(os.path.basename(tsv), sep='\t', index=False)
        EOF

        python3 mergePOLYX.py ~{polyx_vcf} ~{tsv} > stdout
    >>>

    output {
        File filtered_polyx_tsv = basename(tsv)
    }
}