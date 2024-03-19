import numpy as np
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

bed = pd.read_csv(bed_file, sep='\t')
bed = bed[bed.cohort==cohort_prefix].copy()
bed['locus_interval'] = bed.CHROM + ':' + bed.START.astype(str) + '-' + bed.END.astype(str)
bed['interval_size'] = bed.END - bed.START

#split-multi
def split_multi_ssc(mt):
    mt = mt.annotate_rows(num_alleles = mt.alleles.size() ) # Add number of alleles at site before split
    # Now split
    sm = hl.split_multi(mt)
    pl = hl.or_missing(hl.is_defined(sm.PL),
                      (hl.range(0, 3).map(lambda i: hl.min(hl.range(0, hl.len(sm.PL))
       .filter(lambda j: hl.downcode(hl.unphased_diploid_gt_index_call(j), sm.a_index) == hl.unphased_diploid_gt_index_call(i))
       .map(lambda j: sm.PL[j])))))
    split_ds = sm.annotate_entries(GT = hl.downcode(sm.GT, sm.a_index),
                                   AD = hl.or_missing(hl.is_defined(sm.AD), [hl.sum(sm.AD) - sm.AD[sm.a_index], sm.AD[sm.a_index]]),
                                   PL = pl) 
        #GQ = hl.cond(hl.is_defined(pl[0]) & hl.is_defined(pl[1]) & hl.is_defined(pl[2]), hl.gq_from_pl(pl), sm.GQ) )
    mt = split_ds.drop('old_locus', 'old_alleles')
    return mt

test_shard = hl.import_vcf(vep_file, reference_genome='GRCh38')
test_shard = split_multi_ssc(test_shard)
test_shard = test_shard.filter_rows(hl.is_snp(test_shard.alleles[0], test_shard.alleles[1]))

test_shard = test_shard.filter_rows(test_shard.filters.size()==0)
test_shard = test_shard.annotate_entries(AB=test_shard.AD[1]/hl.sum(test_shard.AD))

if test_shard.count_rows()==0:
    pd.DataFrame({'locus': [], 'alleles': [], 'AB': []}).to_csv(output_name, sep='\t', index=False)  
else:
    start_locus = test_shard.head(1).locus.collect()[0]
    end_locus = test_shard.tail(1).locus.collect()[0]
    chrom_str_to_int = {f"{i}": i for i in range(1, 23)} | {'X': 23, 'Y': 24}
    end_chrom = chrom_str_to_int[end_locus.contig.split('chr')[1]]
    end_pos = int(end_locus.position)
    start_chrom = chrom_str_to_int[start_locus.contig.split('chr')[1]]
    start_pos = int(start_locus.position)

    bed['chrom_int'] = bed.CHROM.str.split('chr').str[1].replace({'X': 23, 'Y': 24}).astype(int)
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