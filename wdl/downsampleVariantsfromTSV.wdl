version 1.0

import "removeRegionsfromTSV.wdl" as removeRegions

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
        File remove_regions_bed
        File hg38_reference
        File hg38_reference_dict
        File hg38_reference_fai
        Int chunk_size=100000
        Float snv_scale=1
        Float indel_scale=1
        String jvarkit_docker
        String hail_docker
        String genome_build='GRCh38'
        Boolean prioritize_gnomad=true
        Boolean prioritize_coding=true
        RuntimeAttr? runtime_attr_downsample
    }

    Array[Pair[String, Float]] var_types_scales = zip(['SNV', 'Indel'], [snv_scale, indel_scale])

    scatter (pair in var_types_scales) {
        String var_type = pair.left
        Float scale = pair.right
        call getNumVars {
        input:
            reference_tsv=reference_tsv,
            hail_docker=hail_docker,
            var_type=var_type
        }

        call removeRegions.removeRegionsVariants as removeRegionsVariants {
            input:
            remove_regions_bed=remove_regions_bed,
            full_input_tsv=full_input_tsv,
            hail_docker=hail_docker,
            var_type=var_type,
            genome_build=genome_build,
            runtime_attr_override=runtime_attr_downsample,
            prioritize_coding=prioritize_coding
        }
        
        call downsampleVariantsPython {
            input:
            full_input_tsv=removeRegionsVariants.filtered_tsv,
            hail_docker=hail_docker,
            var_type=var_type,
            chunk_size=chunk_size,
            scale=scale,
            num_variants=getNumVars.num_variants,
            runtime_attr_override=runtime_attr_downsample,
            prioritize_gnomad=prioritize_gnomad,
            prioritize_coding=prioritize_coding
        }

        call convertTSVtoVCF {
            input:
            tsv=downsampleVariantsPython.downsampled_tsv,
            runtime_attr_override=runtime_attr_downsample,
            hail_docker=hail_docker,
            genome_build=genome_build
        }

        call annotatePOLYX {
            input:
            vcf_file=convertTSVtoVCF.output_vcf,
            hg38_reference=hg38_reference,
            hg38_reference_fai=hg38_reference_fai,
            hg38_reference_dict=hg38_reference_dict,
            runtime_attr_override=runtime_attr_downsample,
            jvarkit_docker=jvarkit_docker
        }

        call mergePOLYX {
            input:
            polyx_vcf=annotatePOLYX.polyx_vcf,
            tsv=downsampleVariantsPython.downsampled_tsv,
            runtime_attr_override=runtime_attr_downsample,
            hail_docker=hail_docker
        }
    }

    output {
        File downsampled_tsv_SNV = mergePOLYX.downsampled_polyx_tsv[0]
        File downsampled_tsv_Indel = mergePOLYX.downsampled_polyx_tsv[1]
    }
}

task getNumVars {
    input {
        File reference_tsv
        String hail_docker
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
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        cat <<EOF > get_num_vars.py
        import pandas as pd
        import os
        import sys

        reference_tsv = sys.argv[1]
        var_type = sys.argv[2]

        df = pd.read_csv(reference_tsv, sep='\t')
        df = df[df.TYPE==var_type]
        print(df.shape[0])
        EOF

        python3 get_num_vars.py ~{reference_tsv} ~{var_type} > stdout
    >>>

    output {
        Int num_variants = read_lines(stdout())[0]
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

        df = pd.concat(pd.read_csv(tsv, sep='\t', chunksize=100_000))[['CHROM','POS','REF','ALT']].to_csv('temp.tsv', sep='\t',index=False)
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

task downsampleVariantsPython {
    input {
        File full_input_tsv
        String hail_docker
        String var_type
        Int num_variants
        Int chunk_size
        Float scale
        Boolean prioritize_gnomad
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

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    Int desired_num_variants = ceil(num_variants * scale)
    String file_ext = if sub(basename(full_input_tsv), '\\.gz', '')==basename(full_input_tsv) then '.tsv' else '.tsv.gz'
    String output_name = "~{basename(full_input_tsv, file_ext)}_~{var_type}_downsampled.tsv.gz"

    command <<<
        cat <<EOF > downsample.py
        import pandas as pd
        import numpy as np
        import os
        import sys
        import ast

        full_input_tsv = sys.argv[1]
        var_type = sys.argv[2]
        desired_num_variants = int(sys.argv[3])
        output_name = sys.argv[4]
        chunksize = int(sys.argv[5])
        prioritize_gnomad = ast.literal_eval(sys.argv[6].capitalize())
        prioritize_coding = ast.literal_eval(sys.argv[7].capitalize())

        chunks = []
        for chunk in pd.read_csv(full_input_tsv, sep='\t', chunksize=chunksize):
            chunks.append(chunk)

        df = pd.concat(chunks)
        df = df[df.TYPE==var_type].copy()

        if prioritize_coding:
            keep = df[df.isCoding].copy()
            df = df[~df.isCoding].copy()  # only downsample from noncoding
            desired_num_variants = max(0, desired_num_variants - keep.shape[0])
        else:
            keep = pd.DataFrame()
        if prioritize_gnomad:
            if 'gnomAD_max_AF' not in df.columns:
                df['gnomAD_max_AF'] = df[['gnomADe_AF', 'gnomADg_AF']].max(axis=1)
            in_gnomad = df[~df.gnomAD_max_AF.isna()].reset_index(drop=True)
            not_in_gnomad = df[df.gnomAD_max_AF.isna()].reset_index(drop=True)

            if in_gnomad.shape[0]>desired_num_variants:
                num_per_sample = int(desired_num_variants / in_gnomad.SAMPLE.unique().size)
                in_gnomad = in_gnomad.groupby('SAMPLE').apply(lambda s: s.sample(min(len(s), num_per_sample))).reset_index(drop=True)

            if not not_in_gnomad.empty:
                num_per_sample = int(max(0, desired_num_variants-in_gnomad.shape[0]) / not_in_gnomad.SAMPLE.unique().size)
                not_in_gnomad = not_in_gnomad.groupby('SAMPLE').apply(lambda s: s.sample(min(len(s), num_per_sample))).reset_index(drop=True)
                df = pd.concat([not_in_gnomad, in_gnomad])
            else:
                df = in_gnomad
        else:
            num_per_sample = int(desired_num_variants / df.SAMPLE.unique().size)
            df = df.groupby('SAMPLE').apply(lambda s: s.sample(min(len(s), num_per_sample)))
        df = pd.concat([keep, df])
        df.to_csv(output_name, sep='\t', index=False)
        EOF

        python3 downsample.py ~{full_input_tsv} ~{var_type} ~{desired_num_variants} ~{output_name} ~{chunk_size} \
        ~{prioritize_gnomad} ~{prioritize_coding}
    >>>

    output {
        File downsampled_tsv = output_name
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
        File downsampled_polyx_tsv = basename(tsv)
    }
}