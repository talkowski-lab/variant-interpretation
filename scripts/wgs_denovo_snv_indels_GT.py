#!/usr/bin/env python
# coding: utf-8

import hail as hl
import pandas as pd
import numpy as np
import sys
import os

build = 'GRCh38'
vcf_file = sys.argv[1]
ped_uri = sys.argv[2]
af_threshold = float(sys.argv[3])
header_file = sys.argv[4]
cores = sys.argv[5]
mem = int(np.floor(float(sys.argv[6])))
file_ext = sys.argv[7]

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

prefix = os.path.basename(vcf_file).split(file_ext)[0]

tmp_ped = pd.read_csv(ped_uri, sep='\t')
# check ped number of columns
if len(tmp_ped) > 6:
    tmp_ped = tmp_ped.iloc[:,:6]
tmp_ped.to_csv(f"{prefix}.ped", sep='\t', index=False)

ped_uri = f"{prefix}.ped"
pedigree = hl.Pedigree.read(ped_uri, delimiter='\t')

mt = hl.import_vcf(vcf_file, force_bgz=True, array_elements_required=False, call_fields=[], header_file=header_file, reference_genome='GRCh38')

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

mt = split_multi_ssc(mt)

mt = mt.annotate_rows(info=mt.info.annotate(cohort_AC=mt.info.AC[mt.a_index - 1],
                                           cohort_AF=mt.info.AF[mt.a_index - 1]))

mt = mt.filter_entries(mt.info.cohort_AF <= af_threshold)

tm = hl.trio_matrix(mt, pedigree)

tm_denovo = tm.filter_entries((tm.proband_entry.GT.is_het()) & 
                  (tm.father_entry.GT.is_hom_ref()) & 
                  (tm.mother_entry.GT.is_hom_ref()))

tm_denovo_df = tm_denovo.entries().to_pandas()
tm_denovo_df.to_csv(f"{os.path.basename(vcf_file).split(file_ext)[0]}_denovo_GT_AF_filter.tsv", sep='\t', index=False)