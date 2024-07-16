from pyspark.sql import SparkSession
import hail as hl
import numpy as np
import sys
import ast
import os

vcf_file = sys.argv[1]
prefix = sys.argv[2]
cores = sys.argv[3]  # string
mem = int(np.floor(float(sys.argv[4])))
ped_uri = sys.argv[5]
ac_threshold = int(sys.argv[6])
gnomad_af_threshold = float(sys.argv[7])
am_threshold = float(sys.argv[8])
mpc_threshold = float(sys.argv[9])
gnomad_rec_threshold = float(sys.argv[10])
gnomad_dom_threshold = float(sys.argv[11])
loeuf_v2_threshold = float(sys.argv[12])
loeuf_v4_threshold = float(sys.argv[13])

hl.init()

mt = hl.import_vcf(vcf_file, reference_genome='GRCh38', force_bgz=True, call_fields=[], array_elements_required=False)

header = hl.get_vcf_metadata(vcf_file)
csq_columns = header['info']['CSQ']['Description'].split('Format: ')[1].split('|')

# Output 1: grab ClinVar only
clinvar_mt = mt.filter_rows(hl.set(['Pathogenic', 'Likely_pathogenic']).intersection(hl.set(mt.info.CLNSIG)).size()!=0)

# filter out ClinVar benign
mt = mt.filter_rows((hl.is_missing(mt.info.CLNSIG)) |
    ~(mt.info.CLNSIG[0].matches('Benign') | mt.info.CLNSIG[0].matches('benign')))

# filter PASS
mt = mt.filter_rows(mt.filters.size()==0)

# split VEP CSQ string
mt = mt.annotate_rows(vep=mt.info)
transcript_consequences = mt.vep.CSQ.map(lambda x: x.split('\|'))

transcript_consequences_strs = transcript_consequences.map(lambda x: hl.if_else(hl.len(x)>1, hl.struct(**
                                                       {col: x[i] if col!='Consequence' else x[i].split('&')  
                                                        for i, col in enumerate(csq_columns)}), 
                                                        hl.struct(**{col: hl.missing('str') if col!='Consequence' else hl.array([hl.missing('str')])  
                                                        for i, col in enumerate(csq_columns)})))

mt = mt.annotate_rows(vep=mt.vep.annotate(transcript_consequences=transcript_consequences_strs))
mt = mt.annotate_rows(vep=mt.vep.select('transcript_consequences'))

# TODO: remove
# # filter by AlphaMissense
# mt = mt.annotate_rows(am_pathogenicity=hl.array(hl.set(mt.vep.transcript_consequences.am_pathogenicity.filter(lambda x: x!=''))))
# mt = mt.annotate_rows(am_pathogenicity=hl.or_missing(mt.am_pathogenicity.size()>0, hl.float(mt.am_pathogenicity[0])))
# mt = mt.filter_rows(hl.is_missing(mt.am_pathogenicity) | (mt.am_pathogenicity>=0.56))

# filter out variants containing only these consequences
exclude_csqs = ['intergenic_variant', 'upstream_gene_variant', 'downstream_gene_variant',
                'synonymous_variant', 'coding_sequence_variant', 'sequence_variant']

mt = mt.annotate_rows(all_csqs=hl.set(hl.flatmap(lambda x: x, mt.vep.transcript_consequences.Consequence)),  
                             gnomADg_AF=hl.or_missing(hl.array(hl.set(mt.vep.transcript_consequences.gnomADg_AF))[0]!='', 
                                                  hl.float(hl.array(hl.set(mt.vep.transcript_consequences.gnomADg_AF))[0])),
                             gnomADe_AF=hl.or_missing(hl.array(hl.set(mt.vep.transcript_consequences.gnomADe_AF))[0]!='', 
                                                  hl.float(hl.array(hl.set(mt.vep.transcript_consequences.gnomADe_AF))[0])))
mt = mt.annotate_rows(gnomad_af=hl.max([mt.gnomADg_AF, mt.gnomADe_AF]))
mt = mt.filter_rows(hl.set(exclude_csqs).intersection(mt.all_csqs).size()!=mt.all_csqs.size())

# filter by AC and gnomAD AF
mt = mt.filter_rows(mt.info.cohort_AC<=ac_threshold)
mt = mt.filter_rows((mt.gnomad_af<=gnomad_af_threshold) | (hl.is_missing(mt.gnomad_af)))

# TODO: remove
# # filter splice variants and MODERATE/HIGH impact variants
# splice_vars = ['splice_donor_5th_base_variant', 'splice_region_variant', 'splice_donor_region_variant']
# keep_vars = ['non_coding_transcript_exon_variant']

# mt = mt.filter_rows(hl.any(lambda csq: hl.array(splice_vars + keep_vars).contains(csq), mt.all_csqs) |
#                   (hl.any(lambda impact: hl.array(['HIGH','MODERATE']).contains(impact), 
#                           mt.vep.transcript_consequences.IMPACT)))

def filter_mt(mt):
    '''
    mt: can be trio matrix (tm) or matrix table (mt) but must be transcript-level, not variant-level
    '''
    # filter by Consequence
    exclude_csqs = ['intergenic_variant', 'upstream_gene_variant', 'downstream_gene_variant',
                    'synonymous_variant', 'coding_sequence_variant', 'sequence_variant']
    mt = mt.filter_rows(hl.set(exclude_csqs).intersection(
        hl.set(mt.vep.transcript_consequences.Consequence)).size()==0)

    # filter by Impact and splice/noncoding consequence
    splice_vars = ['splice_donor_5th_base_variant', 'splice_region_variant', 'splice_donor_region_variant']
    keep_vars = ['non_coding_transcript_exon_variant']
    mt = mt.filter_rows(
        (hl.set(splice_vars + keep_vars).intersection(
            hl.set(mt.vep.transcript_consequences.Consequence)).size()>0) |
        (hl.array(['HIGH', 'MODERATE']).contains(
        mt.vep.transcript_consequences.IMPACT))
        )

    # filter by Alpha Missense
    mt = mt.filter_rows(hl.if_else(mt.vep.transcript_consequences.am_pathogenicity=='', 1, 
               hl.float(mt.vep.transcript_consequences.am_pathogenicity))>=am_threshold)
    return mt 

# Phasing
pedigree = hl.Pedigree.read(ped_uri)

# Output 2: OMIM Recessive --> CompHets + Homozygous in probands + XLR in males 
tm = hl.trio_matrix(mt, pedigree, complete_trios=True)
phased_tm = hl.experimental.phase_trio_matrix_by_transmission(tm, call_field='GT', phased_call_field='PBT_GT')

test_phased_tm = phased_tm.filter_entries((hl.is_defined(phased_tm.proband_entry.PBT_GT)) 
                                          & (phased_tm.proband_entry.PBT_GT!=hl.parse_call('0|0')))

gene_phased_tm = test_phased_tm.explode_rows(test_phased_tm.vep.transcript_consequences)
# gene_phased_tm = gene_phased_tm.filter_rows(gene_phased_tm.vep.transcript_consequences.Feature!='')
gene_phased_tm = filter_mt(gene_phased_tm)

# XLR only
xlr_phased_tm = gene_phased_tm.filter_rows(gene_phased_tm.vep.transcript_consequences.OMIM_inheritance_code.matches('4'))

# OMIM recessive only
gene_phased_tm = gene_phased_tm.filter_rows(
    (gene_phased_tm.vep.transcript_consequences.OMIM_inheritance_code.matches('2')) |    # OMIM recessive
    ((hl.is_missing(gene_phased_tm.vep.transcript_consequences.OMIM_inheritance_code)) &  # not OMIM recessive with gnomAD AF and MPC filters
                            ((gene_phased_tm.gnomad_af<=gnomad_rec_threshold) | (hl.is_missing(gene_phased_tm.gnomad_af))) &
                            ((gene_phased_tm.info.MPC>=mpc_threshold) | (hl.is_missing(gene_phased_tm.info.MPC)))
    )
)

potential_comp_hets = (gene_phased_tm.group_rows_by(gene_phased_tm.vep.transcript_consequences.Feature)
    .aggregate_rows(locus_alleles = hl.agg.collect(gene_phased_tm.row_key))
    .aggregate_entries(proband_PBT_GT = hl.agg.collect(gene_phased_tm.proband_entry.PBT_GT))).result()

potential_comp_hets = potential_comp_hets.filter_rows(
    hl.agg.count_where(hl.set(potential_comp_hets.proband_PBT_GT).size()>1)>0
)

potential_comp_hets = potential_comp_hets.explode_rows(potential_comp_hets.locus_alleles)
potential_comp_hets = potential_comp_hets.key_rows_by(potential_comp_hets.locus_alleles.locus, potential_comp_hets.locus_alleles.alleles, 'Feature')
potential_comp_hets = potential_comp_hets.filter_entries(hl.set(potential_comp_hets.proband_PBT_GT).size()>1)

gene_phased_tm = gene_phased_tm.key_rows_by('locus', 'alleles', gene_phased_tm.vep.transcript_consequences.Feature)
gene_phased_tm = gene_phased_tm.annotate_entries(proband_PBT_GT_set=hl.set(
    potential_comp_hets[gene_phased_tm.row_key, gene_phased_tm.col_key].proband_PBT_GT))

gene_phased_tm = gene_phased_tm.filter_entries(gene_phased_tm.proband_PBT_GT_set.size()>1)
omim_rec_comp_hets = gene_phased_tm.semi_join_rows(potential_comp_hets.rows()).key_rows_by('locus', 'alleles')#.count()

omim_rec_hom_var = gene_phased_tm.filter_entries(gene_phased_tm.proband_entry.GT.is_hom_var())
omim_rec_hom_var = omim_rec_hom_var.filter_rows(hl.agg.count_where(
    hl.is_defined(omim_rec_hom_var.proband_entry.GT))>0).key_rows_by('locus', 'alleles')

xlr_phased_tm = xlr_phased_tm.annotate_rows(variant_type='XLR')
omim_rec_hom_var = omim_rec_hom_var.annotate_rows(variant_type='hom_var')
omim_rec_comp_hets = omim_rec_comp_hets.annotate_rows(variant_type='comphet')

omim_rec_merged = xlr_phased_tm.union_rows(omim_rec_hom_var.drop('Feature', 'proband_PBT_GT_set'))\
.union_rows(omim_rec_comp_hets.drop('Feature', 'proband_PBT_GT_set'))
omim_rec_df = omim_rec_merged.entries().to_pandas()

# Output 3: OMIM Dominant
gene_tm = tm.explode_rows(tm.vep.transcript_consequences)
gene_tm = filter_mt(gene_tm)
omim_dom = gene_tm.filter_rows((
        (gene_tm.vep.transcript_consequences.OMIM_inheritance_code.matches('1')) &  # OMIM dominant with gnomAD AF and MPC filters 
            ((gene_tm.gnomad_af<=gnomad_dom_threshold) | (hl.is_missing(gene_tm.gnomad_af))) & 
            ((gene_tm.info.MPC>=mpc_threshold) | (hl.is_missing(gene_tm.info.MPC)))) |
        ((hl.is_missing(gene_tm.vep.transcript_consequences.OMIM_inheritance_code)) &  # not OMIM dominant with LOEUF v2/v4 filters
            (hl.if_else(gene_tm.vep.transcript_consequences.LOEUF_v2=='', 0, 
                hl.float(gene_tm.vep.transcript_consequences.LOEUF_v2))<=loeuf_v2_threshold) | 
            (hl.if_else(gene_tm.vep.transcript_consequences.LOEUF_v4=='', 0, 
                hl.float(gene_tm.vep.transcript_consequences.LOEUF_v4))<=loeuf_v4_threshold)
        )
    )
omim_dom_df = omim_dom.entries().to_pandas()

hl.export_vcf(clinvar_mt, prefix+'_clinvar_variants.vcf.bgz', metadata=header)
omim_rec_df.to_csv(prefix+'_OMIM_recessive.tsv', sep='\t', index=False)
omim_dom_df.to_csv(prefix+'_OMIM_dominant.tsv', sep='\t', index=False)
# hl.export_vcf(splice_mt, prefix+'_splice_variants.vcf.bgz', metadata=header)
# hl.export_vcf(nc_impact_mt, prefix+'_noncoding_high_moderate_impact_variants.vcf.bgz', metadata=header)