version 1.0

import "wes-denovo-helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow getBAFinBED {
    input {
        File bed_file
        Array[File] vep_vcf_files
        String cohort_prefix
        String hail_docker
    }

    scatter (vep_file in vep_vcf_files) {
        call getBAF {
            input:
            bed_file=bed_file,
            vep_file=vep_file,
            cohort_prefix=cohort_prefix,
            hail_docker=hail_docker
        }
    }

    call helpers.mergeResults as mergeResults {
        input:
        tsvs=getBAF.baf_tsv,
        hail_docker=hail_docker,
        merged_filename=cohort_prefix + '_AB_SNVs_locus_intervals.tsv'
    }

    output {
        File merged_baf_tsv = mergeResults.merged_tsv
    }
}

task getBAF {
    input {
        File bed_file
        File vep_file
        String cohort_prefix
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vep_file, "GB")
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

    String output_name = basename(vep_file, '.vcf.bgz') + '_AB_SNVs_locus_intervals.tsv'

    command <<<
    cat <<EOF > get_baf.py
    import pandas as pd
    import os
    import sys
    import hail as hl

    bed_file = sys.argv[1]
    cohort_prefix = sys.argv[2]
    vep_file = sys.argv[3]
    cores = sys.argv[4]
    mem = int(np.floor(float(sys.argv[5])))
    output_name = sys.argv[6]

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{mem}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{mem}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    bed = pd.read_csv(bed_file, sep='\t', header=None)
    bed = bed[bed.cohort==cohort_prefix].copy()
    bed['locus_interval'] = bed.CHROM + ':' + bed.START.astype(str) + '-' + bed.END.astype(str)
    bed['interval_size'] = bed.END - bed.START

    test_shard = hl.import_vcf(vep_file, reference_genome='GRCh38')
    test_shard = split_multi_ssc(test_shard)
    test_shard = test_shard.filter_rows(hl.is_snp(test_shard.alleles[0], test_shard.alleles[1]))

    test_shard = test_shard.filter_rows(test_shard.filters.size()==0)
    test_shard = test_shard.annotate_entries(AB=test_shard.AD[1]/hl.sum(test_shard.AD))

    start_locus = test_shard.head(1).locus.collect()[0]
    end_locus = test_shard.tail(1).locus.collect()[0]
    end_chrom = int(end_locus.contig.split('chr')[1])
    end_pos = int(end_locus.position)
    start_chrom = int(start_locus.contig.split('chr')[1])
    start_pos = int(start_locus.position)

    bed['chrom_int'] = bed.CHROM.str.split('chr').str[1].astype(int)
    bed['not_in_vcf'] = bed.apply(lambda row: (row.chrom_int > end_chrom) | 
                                            ((row.chrom_int == end_chrom) & (row.START > end_pos)) |
                                            (row.chrom_int < start_chrom) | 
                    ((row.chrom_int == start_chrom) & (row.END < start_pos)), axis=1)

    def test_interval(locus_interval, sample, mt):
        test_mt = hl.filter_intervals(mt, [hl.parse_locus_interval(locus_interval, 'GRCh38')])
        print(f"number of SNVs in cohort at {locus_interval}: {test_mt.count_rows()}")

        test_mt = test_mt.filter_cols(test_mt.s==sample)
        test_mt = hl.variant_qc(test_mt)
        test_mt = test_mt.filter_rows(test_mt.variant_qc.AC[1]>0)

        print(f"number of SNVs in sample {sample} at {locus_interval}: {test_mt.count_rows()}")

        test_ab = test_mt.AB.collect()
        test_alleles = test_mt.alleles.collect()
        test_locus = test_mt.locus.collect()
        test_df = pd.DataFrame({'locus': test_locus, 'alleles': test_alleles, 'AB': test_ab})
        test_df['locus_interval'] = locus_interval
        test_df['SAMPLE'] = sample
        return test_df

    test_df = pd.DataFrame()
    for i, row in bed[~bed.not_in_vcf].iterrows():
        test_df = pd.concat([test_df, test_interval(row.locus_interval, row.SAMPLE, test_shard)])  
    test_df.to_csv(output_name, sep='\t', index=False)  
    EOF

    python3 get_baf.py ~{bed_file} ~{cohort_prefix} ~{vep_file} ~{cpu_cores} ~{memory} ~{output_name}
    >>>

    output {
        File baf_tsv = output_name
    }
}