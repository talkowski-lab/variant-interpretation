import pandas as pd
import numpy as np
import hail as hl
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
import ast
from gnomad.resources.grch38 import gnomad
from gnomad.sample_qc import relatedness

vcf_uri = sys.argv[1]
bed_file = sys.argv[2]
cohort_prefix = sys.argv[3]
ped_uri = sys.argv[4]
cores = sys.argv[5]  # string
mem = int(np.floor(float(sys.argv[6])))

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

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

mt = hl.import_vcf(vcf_uri, reference_genome='GRCh38', force_bgz=True, array_elements_required=False)
mt = split_multi_ssc(mt)

# somalier sites
intervals = hl.import_bed(bed_file, reference_genome='GRCh38')
mt = mt.filter_rows(hl.is_defined(intervals[mt.locus]))

rel = hl.pc_relate(mt.GT, 0.01, k=10)

rel = rel.annotate(relationship = relatedness.get_relationship_expr(rel.kin, rel.ibd0, rel.ibd1, rel.ibd2, 
                                                   first_degree_kin_thresholds=(0.19, 0.4), second_degree_min_kin=0.1, 
                                                   ibd0_0_max=0.025, ibd0_25_thresholds=(0.1, 0.425), ibd1_0_thresholds=(-0.15, 0.1), 
                                                   ibd1_50_thresholds=(0.275, 0.75), ibd1_100_min=0.75, ibd2_0_max=0.125, 
                                                   ibd2_25_thresholds=(0.1, 0.5), ibd2_100_thresholds=(0.75, 1.25))
)

rel_df = rel.to_pandas()

ped = pd.read_csv(ped_uri, sep='\t').iloc[:, :6]
ped.columns = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
ped.index = ped.sample_id

vcf_samps = mt.s.collect()

for s, row in ped[ped.role=='Proband'].iterrows():
    proband, father, mother = row[['sample_id', 'paternal_id', 'maternal_id']]
    if father=='0':
        dad_df = rel_df[(rel_df['i.s']==proband) | (rel_df['j.s']==proband)]
#         break
    else:
        ped.loc[s,'father_status'] = rel_df[((rel_df['i.s']==proband) & (rel_df['j.s']==father)) 
           | ((rel_df['i.s']==father) & (rel_df['j.s']==proband))].relationship.item()
    if mother=='0':
        mom_df = rel_df[(rel_df['i.s']==proband) | (rel_df['j.s']==proband)]
#         break
    else:
        ped.loc[s,'mother_status'] = rel_df[((rel_df['i.s']==proband) & (rel_df['j.s']==mother)) 
           | ((rel_df['i.s']==mother) & (rel_df['j.s']==proband))].relationship.item()
