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
ped_uri = sys.argv[7]

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
ped.columns = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
ped.index = ped.sample_id

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

test_shard = hl.filter_intervals(test_shard, [hl.parse_locus_interval(locus_interval, 'GRCh38') for locus_interval in bed.locus_interval])
test_shard = test_shard.filter_rows(hl.is_snp(test_shard.alleles[0], test_shard.alleles[1]))

test_shard = test_shard.filter_rows(test_shard.filters.size()==0)
test_shard = test_shard.annotate_entries(AB=test_shard.AD[1]/hl.sum(test_shard.AD))

roles = ['sample','father','mother']
output_cols = ['locus', 'alleles', 'locus_interval', 'SV_type', 'SAMPLE'] + [f"AB_{role}" for role in roles] + [f"GT_{role}" for role in roles] 

if test_shard.count_rows()==0:
    print('no overlapping SNVs, outputting empty file...')
    pd.DataFrame({col: [] for col in output_cols}).to_csv(output_name, sep='\t', index=False)  
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

    def test_interval(locus_interval, sample, sv_type, mt):
        mt = hl.filter_intervals(mt, [hl.parse_locus_interval(locus_interval, 'GRCh38')])
        print(f"number of SNVs in cohort at {locus_interval}: {mt.count_rows()}")

        father, mother = ped.loc[sample, ['paternal_id','maternal_id']].tolist()
        trio = [sample, father, mother]

        sample_mt = mt.filter_cols(mt.s==sample)
        sample_mt = hl.variant_qc(sample_mt)
        sample_mt = sample_mt.filter_rows(sample_mt.variant_qc.AC[1]>0)

        n_sample_vars = sample_mt.count_rows()
        print(f"number of SNVs in sample {sample} at {locus_interval}: {n_sample_vars}")

        if n_sample_vars==0:
            return pd.DataFrame({col: [] for col in output_cols})
    
        trio_mt = mt.filter_cols((mt.s==sample) | (mt.s==father) | (mt.s==mother))
        trio_mt = trio_mt.filter_rows(hl.is_defined(sample_mt.rows()[trio_mt.row_key]))
        trio_mat = trio_mt.select_rows().select_entries('AB', 'GT').entries().to_pandas()

        sample_roles = {sample: 'sample', father: 'father', mother: 'mother'}

        trio_mat['role'] = trio_mat.s.map(sample_roles)
        trio_mat['alleles'] = trio_mat.alleles.apply(':'.join)

        trio_mat_new = pd.pivot_table(trio_mat, index=['locus','alleles'], columns=['role'], values=['AB','GT'])
        trio_mat_new.columns = trio_mat_new.columns.map('_'.join)

        trio_mat_new = trio_mat_new.reset_index()
        trio_mat_new['locus_interval'] = locus_interval
        trio_mat_new['SAMPLE'] = sample
        trio_mat_new['SV_type'] = sv_type

        return trio_mat_new
    
    test_df = pd.DataFrame()
    for i, row in bed[~bed.not_in_vcf].iterrows():
        test_df = pd.concat([test_df, test_interval(row.locus_interval, row.SAMPLE, row.TYPE, test_shard)])  
    test_df[output_cols].to_csv(output_name, sep='\t', index=False)  