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
window_size = float(sys.argv[8])

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
ped.columns = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
ped.index = ped.sample_id

chroms = [f"chr{x}" for x in range(1, 23)] + ["chrX", "chrY"]
chrom_lengths = {chrom: hl.eval(hl.contig_length(chrom, reference_genome='GRCh38')) for chrom in chroms}

bed = pd.read_csv(bed_file, sep='\t')
bed = bed[bed.cohort==cohort_prefix].copy()
bed['locus_interval'] = bed.CHROM + ':' + bed.START.astype(str) + '-' + bed.END.astype(str)
bed['interval_size'] = bed.END - bed.START
bed['window_start'] = (bed.START - window_size*bed.interval_size).apply(lambda x: max(int(x), 1))
bed['window_end'] = bed.apply(lambda row: min(chrom_lengths[row.CHROM], int(row.END + window_size*row.interval_size)), axis=1)
bed['window_locus_interval'] = bed.CHROM + ':' + bed.window_start.astype(str) + '-' + bed.window_end.astype(str)

#split-multi
def split_multi_ssc(mt):
    mt = mt.annotate_rows(num_alleles = mt.alleles.size() ) # Add number of alleles at site before split
    # only split variants that aren't already split
    bi = mt.filter_rows(hl.len(mt.alleles) == 2)
    bi = bi.annotate_rows(a_index=1, was_split=False, old_locus=bi.locus, old_alleles=bi.alleles)
    multi = mt.filter_rows(hl.len(mt.alleles) > 2)
    # Now split
    split = hl.split_multi(multi, permit_shuffle=True)
    sm = split.union_rows(bi)
    # sm = hl.split_multi(mt, permit_shuffle=True)
    if 'PL' in list(mt.entry.keys()):
        pl = hl.or_missing(hl.is_defined(sm.PL),
                        (hl.range(0, 3).map(lambda i: hl.min(hl.range(0, hl.len(sm.PL))
        .filter(lambda j: hl.downcode(hl.unphased_diploid_gt_index_call(j), sm.a_index) == hl.unphased_diploid_gt_index_call(i))
        .map(lambda j: sm.PL[j])))))
        sm = sm.annotate_entries(PL = pl)
    split_ds = sm.annotate_entries(GT = hl.downcode(sm.GT, sm.a_index),
                                   AD = hl.or_missing(hl.is_defined(sm.AD), [hl.sum(sm.AD) - sm.AD[sm.a_index], sm.AD[sm.a_index]])
                                   ) 
        #GQ = hl.cond(hl.is_defined(pl[0]) & hl.is_defined(pl[1]) & hl.is_defined(pl[2]), hl.gq_from_pl(pl), sm.GQ) )
    mt = split_ds.drop('old_locus', 'old_alleles')
    return mt

mt = hl.import_vcf(vep_file, reference_genome='GRCh38', array_elements_required=False, call_fields=[])
mt = split_multi_ssc(mt)

mt = hl.filter_intervals(mt, [hl.parse_locus_interval(locus_interval, 'GRCh38') for locus_interval in bed.locus_interval])
mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))

mt = mt.filter_rows(mt.filters.size()==0)
mt = mt.annotate_entries(AB=mt.AD[1]/hl.sum(mt.AD))

# filter out regions if sample not in VCF
vcf_samps = mt.s.collect()
bed = bed[bed.SAMPLE.isin(vcf_samps)].copy()

roles = ['sample','father','mother']
output_cols = ['locus', 'alleles', 'window_locus_interval', 'locus_interval', 'SV_type', 'SAMPLE', 'pipeline_id', 'vep_shard'] + [f"AB_{role}" for role in roles] + [f"GT_{role}" for role in roles] 

def test_interval(locus_interval, og_locus_interval, sample, pipeline_id, sv_type, mt):
    mt = hl.filter_intervals(mt, [hl.parse_locus_interval(locus_interval, 'GRCh38')])
    n_cohort_vars = mt.count_rows()
    to_return = {'window_locus_interval': [locus_interval], 'locus_interval': [og_locus_interval], 'SV_type': [sv_type],
                            'SAMPLE': [sample], 'pipeline_id': [pipeline_id],
                            'vep_shard': [os.path.basename(vep_file)]}
    if n_cohort_vars==0:
        print(f"No SNVs in cohort at {locus_interval}...")
        return pd.DataFrame(to_return | {col: [np.nan] for col in np.setdiff1d(output_cols, list(to_return.keys()))})
    print(f"number of SNVs in cohort at {locus_interval}: {n_cohort_vars}")

    father, mother = ped.loc[sample, ['paternal_id','maternal_id']].tolist()

    sample_mt = mt.filter_cols(mt.s==sample)
    sample_mt = hl.variant_qc(sample_mt)
    sample_mt = sample_mt.filter_rows(sample_mt.variant_qc.AC[1]>0)

    n_sample_vars = sample_mt.count_rows()
    if n_sample_vars==0:
        print(f"No SNVs in sample {sample} at {locus_interval}...")
        return pd.DataFrame(to_return | {col: [np.nan] for col in np.setdiff1d(output_cols, list(to_return.keys()))})
    print(f"number of SNVs in sample {sample} at {locus_interval}: {n_sample_vars}")

    trio_mt = mt.filter_cols((mt.s==sample) | (mt.s==father) | (mt.s==mother))
    trio_mt = trio_mt.filter_rows(hl.is_defined(sample_mt.rows()[trio_mt.row_key]))
    trio_mat = trio_mt.select_rows().select_entries('AB', 'GT').entries().to_pandas()

    sample_roles = {sample: 'sample', father: 'father', mother: 'mother'}

    trio_mat['role'] = trio_mat.s.map(sample_roles)
    trio_mat['alleles'] = trio_mat.alleles.apply(':'.join)

    trio_mat_new = trio_mat.pivot(index=['locus','alleles'], columns=['role'], values=['AB','GT'])
    trio_mat_new.columns = trio_mat_new.columns.map('_'.join)

    trio_mat_new = trio_mat_new.reset_index()
    trio_mat_new['window_locus_interval'] = locus_interval
    trio_mat_new['locus_interval'] = og_locus_interval
    trio_mat_new['SAMPLE'] = sample
    trio_mat_new['pipeline_id'] = pipeline_id
    trio_mat_new['SV_type'] = sv_type
    trio_mat_new['vep_shard'] = os.path.basename(vep_file)

    return trio_mat_new

if mt.count_rows()==0:
    print('no overlapping SNVs, outputting empty file...')
    pd.DataFrame({col: [] for col in output_cols}).to_csv(output_name, sep='\t', index=False)  

else:
    start_locus = mt.head(1).locus.collect()[0]
    end_locus = mt.tail(1).locus.collect()[0]
    chrom_str_to_int = {f"{i}": i for i in range(1, 23)} | {'X': 23, 'Y': 24}
    end_chrom = chrom_str_to_int[end_locus.contig.split('chr')[1]]
    end_pos = int(end_locus.position)
    start_chrom = chrom_str_to_int[start_locus.contig.split('chr')[1]]
    start_pos = int(start_locus.position)

    bed['chrom_int'] = bed.CHROM.str.split('chr').str[1].replace({'X': 23, 'Y': 24}).astype(int)
    bed['not_in_vcf'] = bed.apply(lambda row: (row.chrom_int > end_chrom) | 
                                            ((row.chrom_int == end_chrom) & (row.window_start > end_pos)) |
                                            (row.chrom_int < start_chrom) | 
                    ((row.chrom_int == start_chrom) & (row.window_end < start_pos)), axis=1)
    
    test_df = pd.DataFrame()
    tot = bed[~bed.not_in_vcf].shape[0]
    for i, row in bed[~bed.not_in_vcf].reset_index().iterrows():
        if ((i+1)%10==0):
            print(f"Getting SNV BAFs for SV interval {i+1}/{tot}...")
        test_df = pd.concat([test_df, test_interval(row.window_locus_interval, row.locus_interval, row.SAMPLE, row.pipeline_id, row.TYPE, mt)])  
    test_df[np.intersect1d(output_cols, test_df.columns)].to_csv(output_name, sep='\t', index=False)  