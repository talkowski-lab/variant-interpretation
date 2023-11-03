#!/usr/bin/env python
# coding: utf-8

import hail as hl

import pandas as pd
import gcsfs
import os
import sys
from google.cloud import storage

lcr_uri = sys.argv[1]
ped_uri = sys.argv[2]
meta_uri = sys.argv[3]
trio_uri = sys.argv[4]
vcf_uri = sys.argv[5]

build='GRCh38'

mt = hl.import_vcf(vcf_uri, reference_genome=build, force_bgz=True, find_replace=('nul', '.'))

# # Impute Sex

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
    return (mt)

# mt = split_multi_ssc(mt)
# # filter low complexity regions
# lcr = hl.import_bed(lcr_uri, reference_genome=build, force_bgz=True)
# mt = mt.filter_rows(hl.is_defined(lcr[mt.locus]), keep=False)
# # select samples after qc
# meta = hl.import_table(meta_uri).key_by('SampleID')
# mt = mt.annotate_cols(pheno=meta[mt.s])
# mt = mt.filter_cols(mt.pheno.Role != '', keep = True)
# # filter for VQSR - PASS only
# mt = mt.filter_rows(mt.filters.size() == 0, keep = True)
# # filter on depth
# mt = mt.filter_entries( (mt.DP < 10) | (mt.DP > 200), keep = False) 

# # impute sex
# imputed_sex = hl.impute_sex(mt.GT)

filename = os.path.basename(vcf_uri).split('.vcf.gz')[0]
# imputed_sex.flatten().export(f"{filename}_imputed_sex_res.tsv")

# # Relatedness

# ## PC_relate

# kinship = hl.pc_relate(mt.GT, 0.01, k=20, statistics='kin', include_self_kinship=True)

# kinship.flatten().export(f"{filename}_relatedness_pc_relate_res.tsv")


# ## KING
p5k = hl.import_locus_intervals('gs://fc-71e715ea-2fb8-4a20-8560-15ed867dcc7d/resources/regions/reference_files_purcell5k_grch38_liftover_2021-09-14.interval_list_start-1.bed', 
                                 reference_genome='GRCh38') #few variants that are likely most useful (PCA and relatedness)
mt5k = mt.filter_rows(hl.is_defined(p5k[mt.locus]), keep = True)

kinship = hl.king(mt5k.GT)

kinship.flatten().export(f"{filename}_relatedness_king_res.tsv")
