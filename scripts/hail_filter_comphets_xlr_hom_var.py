from pyspark.sql import SparkSession
import hail as hl
import numpy as np
import pandas as pd
import sys
import ast
import os

snv_indel_vcf = sys.argv[1]
clinvar_vcf = sys.argv[2]
sv_vcf = sys.argv[3]
ped_uri = sys.argv[4]
prefix = sys.argv[5]
omim_uri = sys.argv[6]
sv_gene_fields = sys.argv[7].split(',')  # ['PREDICTED_LOF', 'PREDICTED_INTRAGENIC_EXON_DUP']
build = sys.argv[8]
cores = sys.argv[9]  # string
mem = int(np.floor(float(sys.argv[10])))

hl.init(min_block_size=128, 
        local=f"local[*]", 
        spark_conf={
                    "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.speculation": 'true'
                    }, 
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

def filter_mt(mt):
    '''
    mt: can be trio matrix (tm) or matrix table (mt) but must be transcript-level, not variant-level
    '''
    # filter by Consequence
    exclude_csqs = ['intergenic_variant', 'upstream_gene_variant', 'downstream_gene_variant',
                    'synonymous_variant', 'coding_sequence_variant', 'sequence_variant']
    mt = mt.filter_rows(hl.set(exclude_csqs).intersection(
        hl.set(mt.vep.transcript_consequences.Consequence)).size()!=hl.set(mt.vep.transcript_consequences.Consequence).size())

    # filter only canonical transcript or MANE PLUS CLINICAL
    mt = mt.filter_rows((mt.vep.transcript_consequences.CANONICAL=='YES') | 
                        (mt.vep.transcript_consequences.MANE_PLUS_CLINICAL!=''))

    # filter by Impact and splice/noncoding consequence
    splice_vars = ['splice_donor_5th_base_variant', 'splice_region_variant', 'splice_donor_region_variant']
    keep_vars = ['non_coding_transcript_exon_variant']
    mt = mt.filter_rows(
        (hl.set(splice_vars + keep_vars).intersection(
            hl.set(mt.vep.transcript_consequences.Consequence)).size()>0) |
        (hl.array(['HIGH', 'MODERATE']).contains(
        mt.vep.transcript_consequences.IMPACT))
        )
    return mt 

def load_split_vep_consequences(vcf_uri):
    mt = hl.import_vcf(vcf_uri, reference_genome=build, force_bgz=True, call_fields=[], array_elements_required=False)
    csq_columns = hl.get_vcf_metadata(vcf_uri)['info']['CSQ']['Description'].split('Format: ')[1].split('|')

    mt = mt.annotate_rows(vep=mt.info)
    transcript_consequences = mt.vep.CSQ.map(lambda x: x.split('\|'))

    transcript_consequences_strs = transcript_consequences.map(lambda x: hl.if_else(hl.len(x)>1, hl.struct(**
                                                        {col: x[i] if col!='Consequence' else x[i].split('&')  
                                                            for i, col in enumerate(csq_columns)}), 
                                                            hl.struct(**{col: hl.missing('str') if col!='Consequence' else hl.array([hl.missing('str')])  
                                                            for i, col in enumerate(csq_columns)})))

    mt = mt.annotate_rows(vep=mt.vep.annotate(transcript_consequences=transcript_consequences_strs))
    mt = mt.annotate_rows(vep=mt.vep.select('transcript_consequences'))
    return mt

## STEP 1: Merge SNV/Indel VCF with SV VCF (or just one of them)
# Load SNV/Indel VCF
if snv_indel_vcf!='NA':
    snv_mt = load_split_vep_consequences(snv_indel_vcf)

    # Load and merge SNV/Indel ClinVar P/LP VCF
    if clinvar_vcf!='NA':
        clinvar_mt = load_split_vep_consequences(clinvar_vcf) 
        snv_mt = snv_mt.union_rows(clinvar_mt).distinct_by_row()

    # filter SNV/Indel MT
    snv_mt = snv_mt.explode_rows(snv_mt.vep.transcript_consequences)
    snv_mt = filter_mt(snv_mt)

    # filter out empty gene fields
    snv_mt = snv_mt.annotate_rows(gene=snv_mt['vep']['transcript_consequences']['SYMBOL'])
    snv_mt = snv_mt.filter_rows(snv_mt.gene!='')

    snv_mt = snv_mt.annotate_rows(variant_type='SNV/Indel')

# Load SV VCF
if sv_vcf!='NA':
    sv_mt = hl.import_vcf(sv_vcf, reference_genome=build, force_bgz=True, call_fields=[], array_elements_required=False)

    # filter by PASS
    sv_mt = sv_mt.filter_rows((sv_mt.filters.size()==0) | hl.is_missing(sv_mt.filters))

    # ignore genes for CPX SVs
    sv_mt = sv_mt.annotate_rows(gene=hl.or_missing(sv_mt.info.SVTYPE!='CPX', 
                                                   hl.array(hl.set(hl.flatmap(lambda x: x, [sv_mt.info[field] for field in sv_gene_fields])))))
    sv_mt = sv_mt.annotate_rows(variant_type='SV')

    sv_mt = sv_mt.explode_rows(sv_mt.gene)

    # VEP
    if (snv_indel_vcf!='NA'):
        snv_vep_fields = {field: str(snv_mt.vep.transcript_consequences[field].dtype) for field in list(snv_mt.row.vep.transcript_consequences)}
    else:
        snv_vep_fields = {'OMIM_MIM_number': 'array<str>', 'OMIM_inheritance_code': 'str'}
    sv_mt = sv_mt.annotate_rows(vep=hl.struct(transcript_consequences=
            {field: hl.missing(dtype) for field, dtype in snv_vep_fields.items()}))

    # Annotate OMIM in SVs
    omim = hl.import_table(omim_uri).key_by('approvedGeneSymbol')
    sv_mt = sv_mt.key_rows_by('gene')
    sv_mt = sv_mt.annotate_rows(vep=sv_mt.vep.annotate(
        transcript_consequences=sv_mt.vep.transcript_consequences.annotate(
        OMIM_MIM_number=hl.if_else(hl.is_defined(omim[sv_mt.row_key]), omim[sv_mt.row_key].mimNumber, ''),
        OMIM_inheritance_code=hl.if_else(hl.is_defined(omim[sv_mt.row_key]), omim[sv_mt.row_key].inheritance_code, ''))))
    sv_mt = sv_mt.key_rows_by('locus', 'alleles')


if (snv_indel_vcf!='NA') and (sv_vcf!='NA'):
    sv_info_fields, sv_entry_fields = list(sv_mt.row.info), list(sv_mt.entry)
    snv_info_fields, snv_entry_fields = list(snv_mt.row.info), list(snv_mt.entry)

    sv_missing_entry_fields = {field: str(snv_mt[field].dtype) for field in np.setdiff1d(snv_entry_fields, sv_entry_fields)}
    snv_missing_entry_fields = {field: str(sv_mt[field].dtype) for field in np.setdiff1d(sv_entry_fields, snv_entry_fields)}

    sv_missing_info_fields = {field: str(snv_mt.info[field].dtype) for field in np.setdiff1d(snv_info_fields, sv_info_fields)}
    snv_missing_info_fields = {field: str(sv_mt.info[field].dtype) for field in np.setdiff1d(sv_info_fields, snv_info_fields)}

    sv_mt = sv_mt.annotate_entries(**{field: hl.missing(dtype) for field, dtype in sv_missing_entry_fields.items()})
    snv_mt = snv_mt.annotate_entries(**{field: hl.missing(dtype) for field, dtype in snv_missing_entry_fields.items()})

    sv_mt = sv_mt.select_entries(*sorted(list(sv_mt.entry)))
    snv_mt = snv_mt.select_entries(*sorted(list(snv_mt.entry)))

    sv_mt = sv_mt.annotate_rows(info=sv_mt.info.annotate(**{field: hl.missing(dtype) for field, dtype in sv_missing_info_fields.items()}))
    snv_mt = snv_mt.annotate_rows(info=snv_mt.info.annotate(**{field: hl.missing(dtype) for field, dtype in snv_missing_info_fields.items()}))

    sv_mt = sv_mt.annotate_rows(info=sv_mt.info.select(*sorted(list(sv_mt.info))))
    snv_mt = snv_mt.annotate_rows(info=snv_mt.info.select(*sorted(list(snv_mt.info))))

    sv_mt = sv_mt.key_rows_by().select_rows(*sorted(list(sv_mt.row))).key_rows_by('locus','alleles')
    snv_mt = snv_mt.key_rows_by().select_rows(*sorted(list(snv_mt.row))).key_rows_by('locus','alleles')

    # Subset shared samples 
    sv_samps = sv_mt.s.collect()
    snv_samps = snv_mt.s.collect()
    shared_samps = list(np.intersect1d(sv_samps, snv_samps))

    if len(shared_samps)==0:
        shared_samps = ['']

    def align_mt2_cols_to_mt1(mt1, mt2):
        mt1 = mt1.add_col_index()
        mt2 = mt2.add_col_index()
        new_col_order = mt2.index_cols(mt1.col_key).col_idx.collect()
        return mt2.choose_cols(new_col_order)
    
    sv_mt = sv_mt.filter_cols(hl.array(shared_samps).contains(sv_mt.s))
    snv_mt = align_mt2_cols_to_mt1(sv_mt, snv_mt)

    variant_types = 'SV_SNV_Indel'
    merged_mt = sv_mt.union_rows(snv_mt)

elif snv_indel_vcf!='NA':
    variant_types = 'SNV_Indel'
    merged_mt = snv_mt

elif sv_vcf!='NA':
    variant_types = 'SV'
    merged_mt = sv_mt


## STEP 2: Get CompHets
# Mendel errors
def get_mendel_errors(mt, phased_tm):
    all_errors, per_fam, per_sample, per_variant = hl.mendel_errors(mt['GT'], pedigree)
    all_errors_mt = all_errors.key_by().to_matrix_table(row_key=['locus','alleles'], col_key=['s'], row_fields=['fam_id'])
    phased_tm = phased_tm.annotate_entries(mendel_code=all_errors_mt[phased_tm.row_key, phased_tm.col_key].mendel_code)
    return phased_tm

def phase_by_transmission_aggregate_by_gene(tm, mt):
    # filter out calls that are hom ref in proband
    tm = tm.filter_entries(tm.proband_entry.GT.is_non_ref())

    phased_tm = hl.experimental.phase_trio_matrix_by_transmission(tm, call_field='GT', phased_call_field='PBT_GT')
    phased_tm = get_mendel_errors(mt, phased_tm)
    gene_phased_tm = phased_tm.key_rows_by('locus','alleles','gene')

    gene_agg_phased_tm = (gene_phased_tm.group_rows_by(gene_phased_tm.gene)
        .aggregate_rows(locus_alleles = hl.agg.collect(gene_phased_tm.row_key),
                       variant_type = hl.agg.collect(gene_phased_tm.variant_type))
        .aggregate_entries(proband_PBT_GT = hl.agg.collect(gene_phased_tm.proband_entry.PBT_GT).filter(hl.is_defined),
                          proband_GT = hl.agg.collect(gene_phased_tm.proband_entry.GT).filter(hl.is_defined))).result()
    return gene_phased_tm, gene_agg_phased_tm

def get_subset_tm(mt, samples, keep=True, complete_trios=False):
    subset_mt = mt.filter_cols(hl.array(samples).contains(mt.s), keep=keep)

    # remove variants missing in subset samples
    subset_mt = hl.variant_qc(subset_mt)
    subset_mt = subset_mt.filter_rows(subset_mt.variant_qc.AC[1]>0)
    subset_mt = subset_mt.drop('variant_qc')

    subset_tm = hl.trio_matrix(subset_mt, pedigree, complete_trios=complete_trios)
    return subset_tm

def get_non_trio_comphets(mt):
    non_trio_tm = get_subset_tm(mt, trio_samples, keep=False)
    non_trio_gene_phased_tm, non_trio_gene_agg_phased_tm = phase_by_transmission_aggregate_by_gene(non_trio_tm, mt)

    # different criteria for non-trios
    potential_comp_hets_non_trios = non_trio_gene_agg_phased_tm.filter_rows(
            hl.agg.count_where(non_trio_gene_agg_phased_tm.proband_GT.size()>1)>0
    )
          
    potential_comp_hets_non_trios = potential_comp_hets_non_trios.explode_rows(potential_comp_hets_non_trios.locus_alleles)
    potential_comp_hets_non_trios = potential_comp_hets_non_trios.key_rows_by(potential_comp_hets_non_trios.locus_alleles.locus, potential_comp_hets_non_trios.locus_alleles.alleles, 'gene')

    potential_comp_hets_non_trios = potential_comp_hets_non_trios.filter_entries(potential_comp_hets_non_trios.proband_GT.size()>1)
    
    non_trio_gene_phased_tm = non_trio_gene_phased_tm.key_rows_by('locus', 'alleles', 'gene')
    non_trio_gene_phased_tm = non_trio_gene_phased_tm.annotate_entries(proband_GT=
        potential_comp_hets_non_trios[non_trio_gene_phased_tm.row_key, non_trio_gene_phased_tm.col_key].proband_GT,
                                                                      proband_GT_set=hl.set(
        potential_comp_hets_non_trios[non_trio_gene_phased_tm.row_key, non_trio_gene_phased_tm.col_key].proband_GT),
                                                                      proband_PBT_GT_set=hl.set(
        potential_comp_hets_non_trios[non_trio_gene_phased_tm.row_key, non_trio_gene_phased_tm.col_key].proband_PBT_GT))
    non_trio_gene_phased_tm = non_trio_gene_phased_tm.filter_entries(non_trio_gene_phased_tm.proband_GT.size()>1)  # this actually seems necessary but idk why tbh
    gene_phased_tm_comp_hets_non_trios = non_trio_gene_phased_tm.semi_join_rows(potential_comp_hets_non_trios.rows()).key_rows_by('locus', 'alleles')
    return gene_phased_tm_comp_hets_non_trios

def get_trio_comphets(mt):
    trio_tm = get_subset_tm(mt, trio_samples, keep=True, complete_trios=True)

    trio_gene_phased_tm, trio_gene_agg_phased_tm = phase_by_transmission_aggregate_by_gene(trio_tm, mt)

    # different criteria for trios (requires phasing)
    potential_comp_hets_trios = trio_gene_agg_phased_tm.filter_rows(
        hl.agg.count_where(hl.set(trio_gene_agg_phased_tm.proband_PBT_GT).size()>1)>0
    )
    potential_comp_hets_trios = potential_comp_hets_trios.explode_rows(potential_comp_hets_trios.locus_alleles)
    potential_comp_hets_trios = potential_comp_hets_trios.key_rows_by(potential_comp_hets_trios.locus_alleles.locus, potential_comp_hets_trios.locus_alleles.alleles, 'gene')

    potential_comp_hets_trios = potential_comp_hets_trios.filter_entries(hl.set(potential_comp_hets_trios.proband_PBT_GT).size()>1)

    trio_gene_phased_tm = trio_gene_phased_tm.key_rows_by('locus', 'alleles', 'gene')
    trio_gene_phased_tm = trio_gene_phased_tm.annotate_entries(proband_GT=
        potential_comp_hets_trios[trio_gene_phased_tm.row_key, trio_gene_phased_tm.col_key].proband_GT,
                                                               proband_GT_set=hl.set(
        potential_comp_hets_trios[trio_gene_phased_tm.row_key, trio_gene_phased_tm.col_key].proband_GT),
                                                               proband_PBT_GT_set=hl.set(
        potential_comp_hets_trios[trio_gene_phased_tm.row_key, trio_gene_phased_tm.col_key].proband_PBT_GT))

    trio_gene_phased_tm = trio_gene_phased_tm.filter_entries(trio_gene_phased_tm.proband_PBT_GT_set.size()>1)  # this actually seems necessary but idk why tbh
    gene_phased_tm_comp_hets_trios = trio_gene_phased_tm.semi_join_rows(potential_comp_hets_trios.rows()).key_rows_by('locus', 'alleles')
    return gene_phased_tm_comp_hets_trios

def get_transmission(phased_tm_ht):
    phased_tm_ht = phased_tm_ht.annotate(transmission=hl.if_else(phased_tm_ht.proband_entry.PBT_GT==hl.parse_call('0|0'), 'uninherited',
            hl.if_else(phased_tm_ht.proband_entry.PBT_GT==hl.parse_call('0|1'), 'inherited_from_mother',
                        hl.if_else(phased_tm_ht.proband_entry.PBT_GT==hl.parse_call('1|0'), 'inherited_from_father',
                                hl.or_missing(phased_tm_ht.proband_entry.PBT_GT==hl.parse_call('1|1'), 'inherited_from_both'))))
    )
    return phased_tm_ht

# Subset pedigree to samples in VCF, edit parental IDs
vcf_samples = merged_mt.s.collect()
tmp_ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
tmp_ped = tmp_ped[tmp_ped.sample_id.isin(vcf_samples)].copy()
tmp_ped['paternal_id'] = tmp_ped.paternal_id.apply(lambda id: id if id in vcf_samples else '0')
tmp_ped['maternal_id'] = tmp_ped.maternal_id.apply(lambda id: id if id in vcf_samples else '0')
tmp_ped.to_csv(f"{prefix}.ped", sep='\t', index=False)

pedigree = hl.Pedigree.read(f"{prefix}.ped", delimiter='\t')
trio_samples = list(np.intersect1d(vcf_samples,
                              list(np.array([[trio.s, trio.pat_id, trio.mat_id] 
                                             for trio in pedigree.complete_trios() if trio.fam_id!='-9']).flatten())))

## Get CompHets 
# Filter to only in autosomes or PAR
comphet_mt = merged_mt.filter_rows(merged_mt.locus.in_autosome_or_par())
# Filter only OMIM recessive or missing (for SVs)
comphet_mt = comphet_mt.filter_rows((comphet_mt.vep.transcript_consequences.OMIM_inheritance_code.matches('2')) | 
                                    (comphet_mt.vep.transcript_consequences.OMIM_inheritance_code==''))

merged_trio_comphets = get_trio_comphets(comphet_mt)
merged_non_trio_comphets = get_non_trio_comphets(comphet_mt)
merged_trio_comphets = merged_trio_comphets.annotate_cols(trio_status='trio')
merged_non_trio_comphets = merged_non_trio_comphets.annotate_cols(trio_status=
                                                              hl.if_else(merged_non_trio_comphets.fam_id=='-9', 
                                                              'not_in_pedigree', 'non_trio'))
merged_comphets = merged_trio_comphets.entries().union(merged_non_trio_comphets.entries())

# XLR only
merged_tm = hl.trio_matrix(merged_mt, pedigree, complete_trios=False)
gene_phased_tm, gene_agg_phased_tm = phase_by_transmission_aggregate_by_gene(merged_tm, merged_mt)
gene_phased_tm = gene_phased_tm.annotate_cols(trio_status=hl.if_else(gene_phased_tm.fam_id=='-9', 'not_in_pedigree', 
                                                   hl.if_else(hl.array(trio_samples).contains(gene_phased_tm.id), 'trio', 'non_trio')))

xlr_phased_tm = gene_phased_tm.filter_rows((gene_phased_tm.vep.transcript_consequences.OMIM_inheritance_code.matches('4')) |  # OMIM XLR
                                           ((gene_phased_tm.locus.in_x_nonpar()) | (gene_phased_tm.locus.in_x_par())))  # on X chromosome
xlr_phased = xlr_phased_tm.filter_entries((xlr_phased_tm.proband_entry.GT.is_non_ref()) &
                            (~xlr_phased_tm.is_female)).key_rows_by('locus', 'alleles').entries()


# HomVar in proband only
phased_hom_var = gene_phased_tm.filter_entries(gene_phased_tm.proband_entry.GT.is_hom_var())
phased_hom_var = phased_hom_var.filter_entries((phased_hom_var.locus.in_x_nonpar()) &
                            (~phased_hom_var.is_female))  # filter out non-PAR chrX in males
phased_hom_var = phased_hom_var.filter_rows(hl.agg.count_where(
    hl.is_defined(phased_hom_var.proband_entry.GT))>0).key_rows_by('locus', 'alleles').entries()

xlr_phased = xlr_phased.annotate(variant_category='XLR')
phased_hom_var = phased_hom_var.annotate(variant_category='hom_var')
merged_comphets = merged_comphets.annotate(variant_category='comphet')

merged_comphets_xlr_hom_var = merged_comphets.drop('proband_GT','proband_GT_set','proband_PBT_GT_set').union(xlr_phased).union(phased_hom_var)
# Annotate PAR status
merged_comphets_xlr_hom_var = merged_comphets_xlr_hom_var.annotate(in_non_par=~(merged_comphets_xlr_hom_var.locus.in_autosome_or_par()))
merged_comphets_xlr_hom_var = get_transmission(merged_comphets_xlr_hom_var)

merged_comphets_xlr_hom_var.flatten().export(f"{prefix}_{variant_types}_comp_hets_xlr_hom_var.tsv.gz")