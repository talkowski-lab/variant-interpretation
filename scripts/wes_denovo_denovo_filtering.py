import datetime
import pandas as pd
import hail as hl
import numpy as np
import pandas as pd
import ast
import sys

filtered_mt = sys.argv[1]
cohort_prefix = sys.argv[2]
ped_uri = sys.argv[3]
loeuf_file = sys.argv[4]
cores = sys.argv[5]
mem = int(np.floor(float(sys.argv[6])))
bucket_id = sys.argv[7]

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

mt = hl.read_matrix_table(filtered_mt)

# Step 3: De novo filtering

## De novo calling preparation

# Annotate PAR regions (pseudoautosomal regions), hemizygous calls, GQ filter (25), and PL checks

# Add PAR annotations
mt = mt.annotate_rows(XNonPAR = mt.locus.in_x_nonpar(),
                      YNonPAR = mt.locus.in_y_nonpar() )

# Annotate hemizygous calls
mt = mt.annotate_entries(is_hemi = ( (mt.locus.in_x_nonpar() & (mt.is_female == False)) | 
                                      mt.locus.in_y_nonpar() | mt.locus.in_mito() ))

# Apply GQ 25 filter across the board
mt = mt.filter_entries(mt.GQ >= 25, keep = True)

#check for PL missingness
mt = mt.annotate_entries(PL = ( hl.case()
                               .when(mt.GT.is_hom_ref() & hl.is_missing(mt.PL), 
                                    [0, mt.GQ, hl.max(mt.GQ, 3*(mt.AD[0] - mt.AD[1]))] )
                               .default(mt.PL) )
                        )

## De novo calling

# Note: De novo function modified from current Hail source 1) to fix the fact that current code is broken 
# in hemizygous regions and 2) to not assume that child genotype is correct

# Import other things used in https://github.com/hail-is/hail/blob/main/hail/python/hail/methods/family_methods.py
from hail.genetics.pedigree import Pedigree
from hail.matrixtable import MatrixTable
from hail.expr import expr_float64
from hail.table import Table
from hail.typecheck import typecheck, numeric
from hail.methods.misc import require_biallelic

tmp_ped = pd.read_csv(ped_uri, sep='\t')
# check ped number of columns
if len(tmp_ped) > 6:
    tmp_ped = tmp_ped.iloc[:,:6]
    
# subset ped to samples in mt
samps = mt.s.collect()
tmp_ped = tmp_ped[tmp_ped.iloc[:,1].isin(samps)]  # sample_id
tmp_ped = tmp_ped.drop_duplicates('sample_id')    

tmp_ped.to_csv(f"{cohort_prefix}.ped", sep='\t', index=False)

ped_uri = f"{cohort_prefix}.ped"
pedigree = hl.Pedigree.read(ped_uri, delimiter='\t')

#de novo calling script by kyle, version 16

@typecheck(mt=MatrixTable,
           pedigree=Pedigree,
           pop_frequency_prior=expr_float64,
           min_gq=int,
           min_p=numeric,
           max_parent_ab=numeric,
           min_child_ab=numeric,
           min_dp_ratio=numeric,
           ignore_in_sample_allele_frequency=bool)

def kyles_de_novo_v16(mt: MatrixTable,
            pedigree: Pedigree,
            pop_frequency_prior,
            *,
            min_gq: int = 20,
            min_p: float = 0.05,
            max_parent_ab: float = 0.05,
            min_child_ab: float = 0.20,
            min_dp_ratio: float = 0.10,
            ignore_in_sample_allele_frequency: bool = False) -> Table:
  
    MOM_DE_NOVO_PRIOR = 1 / 100000000
    DAD_DE_NOVO_PRIOR = 1 / 25000000
    MITO_DE_NOVO_PRIOR = 10 * MOM_DE_NOVO_PRIOR # Rough estimate, could use review for mitochondrial applications
    
    MIN_POP_PRIOR = 100 / 30000000

    required_entry_fields = {'GT', 'AD', 'DP', 'GQ', 'PL'}
    missing_fields = required_entry_fields - set(mt.entry)
    if missing_fields:
        raise ValueError(f"'de_novo': expected 'MatrixTable' to have at least {required_entry_fields}, "
                         f"missing {missing_fields}")

    pop_frequency_prior = hl.case() \
        .when((pop_frequency_prior >= 0) & (pop_frequency_prior <= 1), pop_frequency_prior) \
        .or_error(hl.str("de_novo: expect 0 <= pop_frequency_prior <= 1, found " + hl.str(pop_frequency_prior)))

    if ignore_in_sample_allele_frequency:
        # this mode is used when families larger than a single trio are observed, in which
        # an allele might be de novo in a parent and transmitted to a child in the dataset.
        # The original model does not handle this case correctly, and so this experimental
        # mode can be used to treat each trio as if it were the only one in the dataset.
        mt = mt.annotate_rows(__prior=pop_frequency_prior,
                              __alt_alleles=hl.int64(1),
                              __site_freq=hl.max(pop_frequency_prior, MIN_POP_PRIOR))
    else:
        # Note that this now requires prior annotation with AC and AN values
        # (did this August 2022 to correct for hemizygous situations)
        n_alt_alleles = mt.AC
        total_alleles = mt.AN
        # subtract 1 from __alt_alleles to correct for the observed genotype
        mt = mt.annotate_rows(__prior=pop_frequency_prior,
                              __alt_alleles=n_alt_alleles,
                              __site_freq=hl.max((n_alt_alleles - 1) / total_alleles,
                                                 pop_frequency_prior,
                                                 MIN_POP_PRIOR))

    mt = require_biallelic(mt, 'de_novo')

    tm = hl.trio_matrix(mt, pedigree, complete_trios=True)
    tm = tm.annotate_rows(__autosomal=tm.locus.in_autosome_or_par(),
                          __x_nonpar=tm.locus.in_x_nonpar(),
                          __y_nonpar=tm.locus.in_y_nonpar(),
                          __mito=tm.locus.in_mito())
    
    autosomal = tm.__autosomal | (tm.__x_nonpar & tm.is_female)
    hemi_x = tm.__x_nonpar & ~tm.is_female
    hemi_y = tm.__y_nonpar & ~tm.is_female
    hemi_mt = tm.__mito

    is_snp = hl.is_snp(tm.alleles[0], tm.alleles[1])
    n_alt_alleles = tm.__alt_alleles
    prior = tm.__site_freq
    
    # Updated for hemizygous child calls to not require a call in uninvolved parent
    has_candidate_gt_configuration = (
        ( autosomal & tm.proband_entry.GT.is_het() & 
          tm.father_entry.GT.is_hom_ref() & tm.mother_entry.GT.is_hom_ref() ) |
        ( hemi_x & tm.proband_entry.GT.is_hom_var() & tm.mother_entry.GT.is_hom_ref() ) |
        ( hemi_y & tm.proband_entry.GT.is_hom_var() & tm.father_entry.GT.is_hom_ref() ) |
        ( hemi_mt & tm.proband_entry.GT.is_hom_var() & tm.mother_entry.GT.is_hom_ref() ) )
    
    failure = hl.missing(hl.tstruct(p_de_novo=hl.tfloat64, confidence=hl.tstr))
    
    kid = tm.proband_entry
    dad = tm.father_entry
    mom = tm.mother_entry

    kid_linear_pl = 10 ** (-kid.PL / 10)
    kid_pp = hl.bind(lambda x: x / hl.sum(x), kid_linear_pl)

    dad_linear_pl = 10 ** (-dad.PL / 10)
    dad_pp = hl.bind(lambda x: x / hl.sum(x), dad_linear_pl)

    mom_linear_pl = 10 ** (-mom.PL / 10)
    mom_pp = hl.bind(lambda x: x / hl.sum(x), mom_linear_pl)

    kid_ad_ratio = kid.AD[1] / hl.sum(kid.AD)
    
    # Try to get these all to an expected value of 0.5
    dp_ratio = (hl.case()
                  .when(hemi_x, kid.DP / mom.DP) # Because mom is diploid but kid is not
                  .when(hemi_y, kid.DP / (2*dad.DP))
                  .when(hemi_mt, kid.DP / (2*mom.DP))
                  .when(tm.__x_nonpar & tm.is_female, (kid.DP / (mom.DP + dad.DP)) * (3/4)) # Because of hemi dad
                  .default( kid.DP / (mom.DP + dad.DP) ))

    def solve(p_de_novo, parent_tot_read_check, parent_alt_read_check):
        return (
            hl.case()
            .when(kid.GQ < min_gq, failure)
            .when((dp_ratio < min_dp_ratio)
                  | (kid_ad_ratio < min_child_ab), failure)
            .when(~parent_tot_read_check, failure)
            .when(~parent_alt_read_check, failure)
            .when(p_de_novo < min_p, failure)
            .when(~is_snp, hl.case()
                  .when((p_de_novo > 0.99) & (kid_ad_ratio > 0.3) & (n_alt_alleles == 1),
                        hl.struct(p_de_novo=p_de_novo, confidence='HIGH'))
                  .when((p_de_novo > 0.5) & (kid_ad_ratio > 0.3) & (n_alt_alleles <= 5),
                        hl.struct(p_de_novo=p_de_novo, confidence='MEDIUM'))
                  .when(kid_ad_ratio > 0.2,
                        hl.struct(p_de_novo=p_de_novo, confidence='LOW'))
                  .or_missing())
            .default(hl.case()
                     .when(((p_de_novo > 0.99) & (kid_ad_ratio > 0.3) & (dp_ratio > 0.2))
                           | ((p_de_novo > 0.99) & (kid_ad_ratio > 0.3) & (n_alt_alleles == 1))
                           | ((p_de_novo > 0.5) & (kid_ad_ratio > 0.3) & (n_alt_alleles < 10) & (kid.DP > 10)),
                           hl.struct(p_de_novo=p_de_novo, confidence='HIGH'))
                     .when((p_de_novo > 0.5) & ((kid_ad_ratio > 0.3) | (n_alt_alleles == 1)),
                           hl.struct(p_de_novo=p_de_novo, confidence='MEDIUM'))
                     .when(kid_ad_ratio > 0.2,
                           hl.struct(p_de_novo=p_de_novo, confidence='LOW'))
                     .or_missing()))
    
    # Call autosomal or pseudoautosomal de novo variants
    def call_auto(kid_pp, dad_pp, mom_pp):
        p_data_given_dn_from_mom = mom_pp[0] * dad_pp[0] * kid_pp[1]
        p_dn_from_mom = ((1 - prior) ** 4) * MOM_DE_NOVO_PRIOR * (1 - DAD_DE_NOVO_PRIOR)
        
        p_data_given_dn_from_dad = mom_pp[0] * dad_pp[0] * kid_pp[1]
        p_dn_from_dad = ((1 - prior) ** 4) * (1 - MOM_DE_NOVO_PRIOR) * DAD_DE_NOVO_PRIOR
        
        p_data_given_mom_alt_allele = mom_pp[1] * dad_pp[0] * kid_pp[1]
        p_mom_alt_allele = prior * ((1 - prior) ** 3) * (1 - MOM_DE_NOVO_PRIOR) * (1 - DAD_DE_NOVO_PRIOR)

        p_data_given_dad_alt_allele = mom_pp[0] * dad_pp[1] * kid_pp[1]
        p_dad_alt_allele = prior * ((1 - prior) ** 3) * (1 - MOM_DE_NOVO_PRIOR) * (1 - DAD_DE_NOVO_PRIOR)

        p_data_given_no_alt_alleles = mom_pp[0] * dad_pp[0] * kid_pp[0]
        p_no_alt_alleles = ((1 - prior) ** 4) * (1 - MOM_DE_NOVO_PRIOR) * (1 - DAD_DE_NOVO_PRIOR)
                               
        p_de_novo = ( (p_data_given_dn_from_mom * p_dn_from_mom) + (p_data_given_dn_from_dad * p_dn_from_dad) ) / (
                      (p_data_given_dn_from_mom * p_dn_from_mom) + (p_data_given_dn_from_dad * p_dn_from_dad) +
                      (p_data_given_mom_alt_allele * p_mom_alt_allele) + (p_data_given_dad_alt_allele * p_dad_alt_allele) +
                      (p_data_given_no_alt_alleles * p_no_alt_alleles) )

        parent_tot_read_check = (hl.sum(mom.AD) > 0) & (hl.sum(dad.AD) > 0)
        parent_alt_read_check = ( ((mom.AD[1] / hl.sum(mom.AD)) <= max_parent_ab) &
                                  ((dad.AD[1] / hl.sum(dad.AD)) <= max_parent_ab) )
            
        return hl.bind(solve, p_de_novo, parent_tot_read_check, parent_alt_read_check)

    # Call hemizygous de novo variants on the X
    def call_hemi_x(kid_pp, mom_pp):
        p_data_given_dn_from_mom = mom_pp[0] * kid_pp[2]
        p_dn_from_mom = ((1 - prior) ** 2) * MOM_DE_NOVO_PRIOR

        p_data_given_mom_alt_allele = mom_pp[1] * kid_pp[2]
        p_mom_alt_allele = prior * (1 - prior) * (1 - MOM_DE_NOVO_PRIOR)

        p_data_given_no_alt_alleles = mom_pp[0] * kid_pp[0]
        p_no_alt_alleles = ((1 - prior) ** 2) * (1 - MOM_DE_NOVO_PRIOR)

        p_de_novo = ( (p_data_given_dn_from_mom * p_dn_from_mom) ) / (
                      (p_data_given_dn_from_mom * p_dn_from_mom) +
                      (p_data_given_mom_alt_allele * p_mom_alt_allele) +
                      (p_data_given_no_alt_alleles * p_no_alt_alleles) )

        parent_tot_read_check = (hl.sum(mom.AD) > 0)
        parent_alt_read_check = (mom.AD[1] / hl.sum(mom.AD)) <= max_parent_ab      

        return hl.bind(solve, p_de_novo, parent_tot_read_check, parent_alt_read_check)

    # Call hemizygous de novo variants on the Y
    def call_hemi_y(kid_pp, dad_pp):
        p_data_given_dn_from_dad = dad_pp[0] * kid_pp[2]
        p_dn_from_dad = (1 - prior) * DAD_DE_NOVO_PRIOR

        p_data_given_dad_alt_allele = dad_pp[2] * kid_pp[2]
        p_dad_alt_allele = prior * (1 - DAD_DE_NOVO_PRIOR)

        p_data_given_no_alt_alleles = dad_pp[0] * kid_pp[0]
        p_no_alt_alleles = (1 - prior) * (1 - DAD_DE_NOVO_PRIOR)
        
        p_de_novo = ( (p_data_given_dn_from_dad * p_dn_from_dad) ) / (
                      (p_data_given_dn_from_dad * p_dn_from_dad) +
                      (p_data_given_dad_alt_allele * p_dad_alt_allele) +
                      (p_data_given_no_alt_alleles * p_no_alt_alleles) )

        parent_tot_read_check = (hl.sum(dad.AD) > 0)
        parent_alt_read_check = (dad.AD[1] / hl.sum(dad.AD)) <= max_parent_ab      

        return hl.bind(solve, p_de_novo, parent_tot_read_check, parent_alt_read_check)       
    
    # Call de novo variants on the mitochondrial chromosome 
    # (treating it as a hemizygous situation and ignoring mito heteroplasmy)
    def call_hemi_mt(kid_pp, mom_pp):
        p_data_given_dn_from_mom = mom_pp[0] * kid_pp[2]
        p_dn_from_mom = (1 - prior) * MITO_DE_NOVO_PRIOR

        p_data_given_mom_alt_allele = mom_pp[2] * kid_pp[2]
        p_mom_alt_allele = prior * (1 - MITO_DE_NOVO_PRIOR)

        p_data_given_no_alt_alleles = mom_pp[0] * kid_pp[0]
        p_no_alt_alleles = (1 - prior) * (1 - MITO_DE_NOVO_PRIOR)
        
        p_de_novo = ( (p_data_given_dn_from_mom * p_dn_from_mom) ) / (
                      (p_data_given_dn_from_mom * p_dn_from_mom) +
                      (p_data_given_mom_alt_allele * p_mom_alt_allele) +
                      (p_data_given_no_alt_alleles * p_no_alt_alleles) )

        parent_tot_read_check = (hl.sum(mom.AD) > 0)
        parent_alt_read_check = (mom.AD[1] / hl.sum(mom.AD)) <= max_parent_ab      

        return hl.bind(solve, p_de_novo, parent_tot_read_check, parent_alt_read_check)       

    
    # Main routine
    de_novo_call = (
        hl.case()
        .when(~has_candidate_gt_configuration, failure)
        .when(autosomal, hl.bind(call_auto, kid_pp, dad_pp, mom_pp))
        .when(hemi_x, hl.bind(call_hemi_x, kid_pp, mom_pp))
        .when(hemi_y, hl.bind(call_hemi_y, kid_pp, dad_pp))
        .when(hemi_mt, hl.bind(call_hemi_mt, kid_pp, mom_pp))
        .or_missing())

    tm = tm.annotate_entries(__call=de_novo_call)
    tm = tm.filter_entries(hl.is_defined(tm.__call))
    entries = tm.entries()
    
    return (entries.select('__site_freq',
                           'proband',
                           'father',
                           'mother',
                           'proband_entry',
                           'father_entry',
                           'mother_entry',
                           'is_female',
                           **entries.__call)
            .rename({'__site_freq': 'prior'}))

# Call de novos using my custom function
de_novo_results = kyles_de_novo_v16(mt, pedigree, pop_frequency_prior = mt.gnomad_non_neuro_AF,
                             max_parent_ab = 0.05, min_child_ab = 0.25, min_dp_ratio = 0.1, min_gq = 25)

# TODO: output (all below)
de_novo_results.flatten().export(f"{cohort_prefix}_wes_final_denovo.txt")

if bucket_id == 'false':
    de_novo_results.write(f"{cohort_prefix}_wes_final_denovo.ht", overwrite = True)
else:
    mt_uris = []
    filename = f"{bucket_id}/hail/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H:%M'))}/{cohort_prefix}_wes_final_denovo.ht"
    mt_uris.append(filename)
    de_novo_results.write(filename, overwrite=True)

# LOEUF annotations

mt = mt.semi_join_rows(de_novo_results.key_by('locus', 'alleles'))
mt.rows().flatten().export(f"{cohort_prefix}_wes_final_denovo_vep.txt")

vep_res = pd.read_csv(f"{cohort_prefix}_wes_final_denovo_vep.txt", sep='\t')
vep_res['alleles'] = vep_res.alleles.apply(ast.literal_eval)
vep_res['ID'] = vep_res.locus + ':' + vep_res.alleles.str.join(':')
vep_res.columns = vep_res.columns.str.replace('info.', '')

loeuf = pd.read_csv(loeuf_file, sep='\t')
loeuf.index = loeuf.gene_name

def get_genes_csq(csq):
    genes = []
    for ind_csq in csq:
        gene = ind_csq.split('|')[3]
        if gene!='':
            genes.append(gene)
    return list(set(genes))

vep_res['CSQ'] = vep_res.CSQ.replace({np.nan: "[]"}).apply(ast.literal_eval)
vep_res['all_genes'] = vep_res.CSQ.apply(get_genes_csq)

all_genes = vep_res.all_genes.apply(pd.Series).stack().unique()

loeuf_vals = loeuf.loc[np.intersect1d(loeuf.index, all_genes), 'LOEUF'].to_dict()
loeuf_tile_vals = loeuf.loc[np.intersect1d(loeuf.index, all_genes), 'LOEUF_tile'].to_dict()

vep_res['LOEUF'] = vep_res.all_genes.apply(lambda gene_list: pd.Series(gene_list).map(loeuf_vals).min())
vep_res['LOEUF_tile'] = vep_res.all_genes.apply(lambda gene_list: pd.Series(gene_list).map(loeuf_tile_vals).min())

vep_res.to_csv(f"{cohort_prefix}_wes_final_denovo_vep.txt", sep='\t', index=False)

## TDT

# Run Hail's built-in TDT function
tdt_table = hl.transmission_disequilibrium_test(mt, pedigree)

# Add TDT results as row annotations
mt = mt.annotate_rows(t = tdt_table.index(mt.row_key).t,
                      u = tdt_table.index(mt.row_key).u)

# Turn it into a trio dataset
td = hl.trio_matrix(mt, pedigree, complete_trios = True)

# TODO: output?
# Write trio dataset
if bucket_id == 'false':
    td.write(f"{cohort_prefix}_trio_tdt.mt", overwrite=True)
else:
    filename = f"{bucket_id}/hail/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H:%M'))}/{cohort_prefix}_trio_tdt.mt"
    mt_uris.append(filename)
    td.write(filename, overwrite=True)

# Parent-aware TDT annotations, plus extras for use in genotype counts
#
# This draws heavily from Jack Kosmicki's function, 
# https://hail.is/docs/0.2/methods/genetics.html#hail.methods.transmission_disequilibrium_test
#
# Requires trio dataset and, like Hail's TDT function, does not cover the Y
#
# Also note:
# 1) Uses sexes from ped and assumes fathers are male and mothers are female
# 2) Does not allow fathers to be het when they should be hemizygous
# 3) To match Jack's function, requires fathers to have a genotype even when considering regions where 
#    proband is hemizygous
def parent_aware_t_u_annotations_v3(td):

    # First decide copy state
    td = td.annotate_entries(autosomal_copy_state = ( hl.case()
            .when(td.locus.in_autosome() | td.locus.in_x_par() | (td.is_female == True), True)
            .when(td.locus.in_x_nonpar() & (td.is_female == False), False)
            .default(hl.missing('bool')) ) )
    # Note: the above uses the "is_female" from the ped and not from the dataset itself
    
    # Now annotate t & u values
    td = td.annotate_entries(
        
        t_from_dad = hl.if_else(
           ( (td.autosomal_copy_state == True) & td.father_entry.GT.is_het() & (td.locus.in_x_nonpar() == False) & 
            ( (td.proband_entry.GT.is_het() & td.mother_entry.GT.is_hom_ref()) | 
              (td.proband_entry.GT.is_hom_var() & 
               (td.mother_entry.GT.is_het() | td.mother_entry.GT.is_hom_var()) )) ), 1, 0),
        
        t_from_mom = hl.if_else(
           ( (td.autosomal_copy_state == True) & td.mother_entry.GT.is_het() &
            ( (td.proband_entry.GT.is_het() & td.father_entry.GT.is_hom_ref()) | 
              (td.proband_entry.GT.is_hom_var() & 
               ((td.father_entry.GT.is_het() & (td.locus.in_x_nonpar() == False)) | td.father_entry.GT.is_hom_var()) )) ) |
           ( (td.autosomal_copy_state == False) & td.mother_entry.GT.is_het() &
               td.proband_entry.GT.is_hom_var() & 
               (td.father_entry.GT.is_hom_ref() | td.father_entry.GT.is_hom_var()) ), 1, 0),
        # I could consider removing any reference at all to father's genotype in this last line    
        
        u_from_dad = hl.if_else(
           ( (td.autosomal_copy_state == True) & td.father_entry.GT.is_het() & (td.locus.in_x_nonpar() == False) & 
            ( (td.proband_entry.GT.is_het() & td.mother_entry.GT.is_hom_var()) | 
              (td.proband_entry.GT.is_hom_ref() & 
               (td.mother_entry.GT.is_het() | td.mother_entry.GT.is_hom_ref()) )) ), 1, 0),
        
        u_from_mom = hl.if_else(
           ( (td.autosomal_copy_state == True) & td.mother_entry.GT.is_het() &
            ( (td.proband_entry.GT.is_het() & td.father_entry.GT.is_hom_var()) | 
              (td.proband_entry.GT.is_hom_ref() & 
               ((td.father_entry.GT.is_het() & (td.locus.in_x_nonpar() == False)) | td.father_entry.GT.is_hom_ref()) )) ) |
           ( (td.autosomal_copy_state == False) & td.mother_entry.GT.is_het() &
               td.proband_entry.GT.is_hom_ref() & 
               (td.father_entry.GT.is_hom_ref() | td.father_entry.GT.is_hom_var()) ), 1, 0),   
        # Again, could consider removing any reference at all to father's genotype in this last line
        
        t_indeterminate = hl.if_else( 
            (td.autosomal_copy_state == True) & td.proband_entry.GT.is_het() & (td.locus.in_x_nonpar() == False) & 
             td.father_entry.GT.is_het() & td.mother_entry.GT.is_het(), 1, 0),

        u_indeterminate = hl.if_else( 
            (td.autosomal_copy_state == True) & td.proband_entry.GT.is_het() & (td.locus.in_x_nonpar() == False) & 
             td.father_entry.GT.is_het() & td.mother_entry.GT.is_het(), 1, 0)
    )
    return (td)

# Add the parent-aware annotations from above
td = parent_aware_t_u_annotations_v3(td)

# TODO: output?
# Write these results to check against built-in TDT function
if bucket_id == 'false':
    td.write(f"{cohort_prefix}_parent_aware_trio_tdt.mt", overwrite = True)
else:
    filename = f"{bucket_id}/hail/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H:%M'))}/{cohort_prefix}_parent_aware_trio_tdt.mt"
    mt_uris.append(filename)
    pd.Series(mt_uris).to_csv('mt_uri.txt', index=False, header=None)
    td.write(filename, overwrite=True)


# TODO: idk the point of the below code...
# # Apply frequency thresholds
# # Apply 0.1% frequency filter on non-neuro gnomad
# td = td.filter_rows(hl.is_missing(td.gnomad_non_neuro_AF) | 
#                    (td.gnomad_non_neuro_AF <= 0.001), keep = True)

# # Apply 0.3% frequency filter on this data
# td = td.filter_rows((td.AC / td.AN) <= 0.003, keep = True)

# # Save rows to export results from Hail's TDT function
# td.rows().write('gs://nd_qc/wes_qc/data_mts/tdt_temp1.ht')

# # Save entries to export results from my own TDT function
# td.filter_entries( (td.t_from_dad + td.t_from_mom + td.t_indeterminate +
#                     td.u_from_dad + td.u_from_mom + td.u_indeterminate) > 0, 
#                   keep = True).entries().write('gs://nd_qc/wes_qc/data_mts/tdt_temp2.ht', overwrite=True)

# # Export row table
# ht = hl.read_table('gs://nd_qc/wes_qc/data_mts/tdt_temp1.ht')

# ht = ht.drop('filters', 'YNonPAR')

# ht.flatten().export('gs://nd_qc/wes_qc/data_mts/tdt_rows.txt.bgz')

# # Export entry table
# ht = hl.read_table('gs://nd_qc/wes_qc/data_mts/tdt_temp2.ht')

# ht = ht.annotate(proband_entry = ht.proband_entry.drop('SB'),
#                  father_entry = ht.father_entry.drop('SB'),
#                  mother_entry = ht.mother_entry.drop('SB', 'is_hemi'))

# ht = ht.drop('filters', 'isSNP', 'YNonPAR')

# ht.flatten().export('gs://nd_qc/wes_qc/data_mts/tdt_entries.txt.bgz')