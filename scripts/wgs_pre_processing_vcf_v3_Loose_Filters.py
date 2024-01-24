#!/usr/bin/env python
# coding: utf-8

# Description: Select only mendelian errors and other basic filters using hail 0.2; Aim to reduce vcf size before splitting by trio unit.
# Shan Dong 2021-5-1
# Update 2021-6-7 (version3): Adding HQ filters (Loose filters); cohortAC filter set to 50; 
# Update 2021-8-18 (version3): cohortAC filter set to 20; Fix sample check step (swtich "sample2" column to "Role" column -- more suitable for the new data sample file)
"""
Usage on googld dataproc: 
hailctl dataproc submit your-dataproc wgs_pre_processing_vcf.py --zone=us-east4-b
"""
import hail as hl

import pandas as pd
import os
import sys
import numpy as np 

build = "GRCh38"
lcr_uri = sys.argv[1]
ped_uri = sys.argv[2]
meta_uri = sys.argv[3]
trio_uri = sys.argv[4]
vcf_uri = sys.argv[5]
header_file = sys.argv[6]
cores = sys.argv[7]  # string
mem = int(np.floor(float(sys.argv[8])))

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

def trim_vcf(vcf_uri, lcr_uri, ped_uri, meta_uri, trio_uri, header_file, vcf_out_uri, build):
    """trim vcf by LCR; select samples after qc (optional); select ME with code = 2 (Fa:HomRef;Mo:HomRef;Child:Het) only; remove loci with low depth"""
    # load vcf
    mt = hl.import_vcf(vcf_uri, array_elements_required=False, reference_genome=build, force_bgz=True, call_fields=[], header_file=header_file, find_replace=('nul', '.'))
    mt = split_multi_ssc(mt)
    # annotate cohort ac to INFO field; cohortAC filter set to 100 
    mt = mt.annotate_rows(info=mt.info.annotate(cohort_AC=mt.info.AC[mt.a_index - 1]))
    mt = mt.filter_rows(mt.info.cohort_AC < 20 , keep=True)
    # filter low complexity regions
    try:
        lcr = hl.import_bed(lcr_uri, reference_genome=build)
    except Exception as e:
        lcr = hl.import_bed(lcr_uri, reference_genome=build, force_bgz=True)
    mt = mt.filter_rows(hl.is_defined(lcr[mt.locus]), keep=False)
    # select samples after qc
    meta = hl.import_table(meta_uri, types={'SampleID': hl.tstr, 'FamID': hl.tstr, 'Role': hl.tstr})
    meta = meta.annotate(s=meta.SampleID).key_by('s')
    mt = mt.annotate_cols(pheno=meta[mt.s])
    mt = mt.filter_cols(mt.pheno.Role != '', keep = True)
    # only keep mendelian errors 
    # mendelian errors code == 2 get parents hom_ref and children het loci
    pedigree = hl.Pedigree.read(ped_uri, delimiter='\t')
    all_errors, per_fam, per_sample, per_variant = hl.mendel_errors(mt['GT'], pedigree)
    mt = mt.semi_join_rows(all_errors.filter(all_errors.mendel_code == 2).key_by('locus', 'alleles'))
    # select specific entries for each variant
    trios = hl.import_table(trio_uri, key='SampleID')
    all_errors = all_errors.drop('fam_id', 'mendel_code')
    all_errors_with_fam = all_errors.key_by('s').join(trios.key_by('SampleID')).drop('FamID', 'TrioID').key_by('locus','alleles')
    all_errors_with_fam = all_errors_with_fam.select(Samples_to_keep = hl.set([all_errors_with_fam.s,all_errors_with_fam.FatherID,all_errors_with_fam.MotherID]))
    all_errors_with_fam = all_errors_with_fam.collect_by_key()
    all_errors_with_fam = all_errors_with_fam.select(Samples_to_keep_all = hl.flatten(hl.set(all_errors_with_fam.values.Samples_to_keep)))
    mt = mt.filter_entries(all_errors_with_fam[mt.row_key].Samples_to_keep_all.contains(mt.s),keep=True)
    # # filter for VQSR - PASS for SNVs only
    mt = mt.filter_rows((hl.is_snp(mt.alleles[0], mt.alleles[1]) &(mt.filters.size() == 0))
                       | hl.is_indel(mt.alleles[0], mt.alleles[1]))
    # filter on depth
    mt = mt.filter_entries( (mt.DP < 10) | (mt.DP > 200), keep = False) 
    # row (variant INFO) level filters - GATK recommendations for short variants
    dn_snv_cond_row = (hl.is_snp(mt.alleles[0], mt.alleles[1])
                        & (mt.qual >= 150)
                        & (mt.info.SOR <= 2.5)
                        & (mt.info.ReadPosRankSum >= -1.4)
                        & (mt.info.QD >= 3.0)
                        & (mt.info.MQ >= 50))
    dn_indel_cond_row = (hl.is_indel(mt.alleles[0], mt.alleles[1])
                          & (mt.qual >= 150)
                          & (mt.info.SOR <= 3)
                          & (mt.info.ReadPosRankSum >= -1.7)
                          & (mt.info.QD >= 4.0)
                          & (mt.info.MQ >= 50))
    mt = mt.filter_rows(dn_snv_cond_row | dn_indel_cond_row, keep = True)
    # parents filters - homozygous
    ab = mt.AD[1]/hl.sum(mt.AD)
    hom_snv_parents = (hl.is_snp(mt.alleles[0], mt.alleles[1])
                    & mt.GT.is_hom_ref()
                    & (mt.GQ >= 30.0)
                    & (ab <= 0.05))
    hom_indel_parents = (hl.is_indel(mt.alleles[0], mt.alleles[1])
                      & mt.GT.is_hom_ref()
                      & (mt.DP >= 16)
                      & (mt.GQ >= 30.0)
                      & (ab <= 0.05))
    # child filters - heterzygous
    het_snv_cond = (hl.is_snp(mt.alleles[0], mt.alleles[1])
                    & mt.GT.is_het()
                    & (mt.GQ >= 99.0)
                    & (ab >= 0.22)
                    & (ab <= 0.78))
    het_indel_cond = (hl.is_indel(mt.alleles[0], mt.alleles[1])
                      & mt.GT.is_het()
                      & (mt.GQ >= 99.0)
                      & (ab >= 0.20)
                      & (ab <= 0.80))
    filter_condition = (hom_snv_parents | hom_indel_parents | het_snv_cond | het_indel_cond)
    mt = mt.filter_entries(filter_condition, keep = True)
    # GQ mean filters
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.gq_stats.mean >= 50, keep = True)
    # clean-up: remove AC = 0 loci
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.AC[1] > 0, keep = True)
    # write to output vcf
    mt = mt.drop('variant_qc')
    hl.export_vcf(mt, vcf_out_uri)
    # header = hl.get_vcf_metadata(vcf_uri) 
    # hl.export_vcf(mt, vcf_out_uri, metadata=header)

vcf_out_uri = os.path.basename(vcf_uri).split('.vcf.gz')[0] + '.preprocessed.vcf.bgz'
trim_vcf(vcf_uri, lcr_uri, ped_uri, meta_uri, trio_uri, header_file, vcf_out_uri, build)

