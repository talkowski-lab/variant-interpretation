import hail as hl
import numpy as np
import sys

vcf_file = sys.argv[1]
cohort_prefix = sys.argv[2]
bucket_id = sys.argv[3]
ped_uri = sys.argv[4]
gnomad_ht_uri = sys.argv[5]
mpc_dir = sys.argv[6]
mpc_chr22_file = sys.argv[7]
purcell5k = sys.argv[8]
cores = sys.argv[9]
mem = int(np.floor(float(sys.argv[10])))

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

mt = hl.import_vcf(vcf_file, reference_genome = 'GRCh38', force_bgz=True)

# Step 1: Annotations

## Filter out discrepant samples from Somalier
# ped = hl.import_table(ped_uri, impute=True).key_by('sample_id')
# mt = mt.filter_cols(hl.is_defined(ped[mt.s]))

## gnomAD exome frequency annotations
gnomad_ht = hl.read_table(gnomad_ht_uri)

# annotate with gnomAD exome frequencies; use "non-neuro" allele frequencies 
mt = mt.annotate_rows(gnomad_non_neuro_AF = 
                      gnomad_ht.index(mt.row_key).freq[hl.eval(gnomad_ht.freq_index_dict["non_neuro"])].AF)

## MPC annotations
mpc = hl.read_table(mpc_dir)

mpc = mpc.annotate(allele_array = [mpc.alleles[0], mpc.alleles[1]])
mpc = mpc.annotate(locus = mpc.locus_38)
mpc = mpc.key_by(mpc.locus, mpc.allele_array)

mpc_chr22 = hl.import_table(mpc_chr22_file, types={"locus": hl.tlocus("GRCh38"), "alleles": hl.tarray(hl.tstr)})
mpc_chr22 = mpc_chr22.annotate(allele_array = [mpc_chr22.alleles[0], mpc_chr22.alleles[1]])
mpc_chr22 = mpc_chr22.key_by(mpc_chr22.locus, mpc_chr22.allele_array)
mpc_chr22 = mpc_chr22.annotate(MPC = hl.float64(mpc_chr22.MPC))

merged_mpc = mpc.select(mpc.MPC).union(mpc_chr22.select(mpc_chr22.MPC))

mt = mt.annotate_rows(MPC=merged_mpc[mt.locus, mt.alleles].MPC)

## pAB annotations

# Add pAB entry field for downstream filtering -- note, only for hets
mt = mt.annotate_entries(pAB = hl.or_missing(mt.GT.is_het(), 
                                             hl.binom_test(mt.AD[1], hl.sum(mt.AD), 0.5, 'two-sided')))

# Run and export Hail's built-in sample QC function
hl.sample_qc(mt, 'sample_qc').cols().flatten().export(
    f"{cohort_prefix}_wes_post_annot_sample_QC_info.txt")

## PCA on 5k variants only
# Create filtered mt
# Load lifted-over Purcell 5k interval list
p5k = hl.import_locus_intervals(purcell5k, 
                                 reference_genome='GRCh38') #few variants that are likely most useful (PCA and relatedness)
mt5k = mt.filter_rows(hl.is_defined(p5k[mt.locus]), keep = True)

#5k variants only
eigenvalues_5k, score_table_5k, loading_table_5k = hl.hwe_normalized_pca(mt5k.GT, k=20, compute_loadings=True)

## also export unfiltered PCA results as mts
score_table_5k.write(f"{cohort_prefix}_wes_pca_score_table_5k.ht", overwrite = True)
loading_table_5k.write(f"{cohort_prefix}_wes_pca_loading_table_5k.ht", overwrite = True)


# Step 2: Filtering

## Basic genotype filtering

# Following filters:

# - All calls: Depth below 7 or above 1000
# - HomRef calls: GQ below 25
# - HomVar calls: PL(HomRef) below 25
# - Het calls:PL(HomRef) below 25

mt = mt.filter_rows((hl.len(mt.filters) == 0))

mt = mt.filter_entries( 
        (mt.DP < 7) | 
        (mt.DP > 1000) |
        ((mt.GT.is_hom_ref()) & (mt.GQ < 25) ) |
        ((mt.GT.is_hom_var()) & (mt.PL[0] < 25) ) |
        ((mt.GT.is_het()) & (mt.PL[0] < 25) ), 
          keep = False
        )

## Remove sites with no variant samples

mt = hl.variant_qc(mt)

mt = mt.filter_rows(mt.variant_qc.AC[1] > 0, keep = True) #applying the above filter to remove those entries where all variants are removed

mt = mt.drop('variant_qc')

## Sex-specific genotype filtering

# Filter:
# - Any call in a female with depth less than 10
# - Any call on the Y in females
# - Any autosomal or PAR call in a male with depth less than 10
# - Het calls in males in hemizygous regions

ped = hl.import_table(ped_uri, impute=True)

original_cols = list(ped.row.keys())
new_cols = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']

ped = ped.rename({old: new for old, new in zip(original_cols, new_cols)})

ped = ped.key_by('sample_id')

mt = mt.annotate_cols(reported_sex = ped[mt.s].sex)

mt = mt.annotate_cols(is_female = (mt.reported_sex == 2))

mt = mt.filter_entries( 
        ((mt.is_female == True) & (mt.DP < 10) ) |
        ((mt.is_female == True) & (mt.locus.in_y_nonpar()) ) |
        ((mt.is_female == False) & (mt.locus.in_autosome_or_par()) & (mt.DP < 10) ) |
        ((mt.is_female == False) & (mt.GT.is_het()) & ( mt.locus.in_x_nonpar() | mt.locus.in_y_nonpar() ) ), 
           keep = False)

# Remove any sites with no variant samples
mt = hl.variant_qc(mt)
mt = mt.filter_rows(mt.variant_qc.AC[1] > 0, keep = True)
mt = mt.drop('variant_qc') # need to recalculate later

## Miscellaneous filters

# filter for het allele balance of 0.25 or more
mt = mt.filter_entries(mt.GT.is_het() & (mt.AD[1] < (0.25 * mt.DP)), keep = False)
mt = hl.variant_qc(mt)
mt = mt.filter_rows(mt.variant_qc.AC[1] > 0, keep = True)
mt = mt.drop('variant_qc')

# filter Het pAB of 0
mt = mt.filter_entries(mt.GT.is_het() & (mt.pAB < 0.000000001), keep = False)
mt = hl.variant_qc(mt)
mt = mt.filter_rows(mt.variant_qc.AC[1] > 0, keep = True)
mt = mt.drop('variant_qc')

# keep only Het informative reads
mt = mt.filter_entries(mt.GT.is_het() & (hl.sum(mt.AD) < (0.9 * mt.DP)), keep = False)
mt = hl.variant_qc(mt)
mt = mt.filter_rows(mt.variant_qc.AC[1] > 0, keep = True)
mt = mt.drop('variant_qc')

# keep only Homvar informative reads
mt = mt.filter_entries(mt.GT.is_hom_var() & (mt.AD[1] < (0.9 * mt.DP)), keep = False)
mt = hl.variant_qc(mt)
mt = mt.filter_rows(mt.variant_qc.AC[1] > 0, keep = True)
mt = mt.drop('variant_qc')

## Sex-aware variant calls

# Define sex-aware variant call rate calculation with Hardy-Weinberg p value

def sex_aware_variant_annotations_with_pHWE(mt):
    num_males = mt.aggregate_cols(hl.agg.count_where(mt.is_female == False))
    num_females = mt.aggregate_cols(hl.agg.count_where(mt.is_female == True))
    
    mt = mt.annotate_rows(
        male_hets = hl.agg.count_where(mt.GT.is_het() & (mt.is_female == False)),
        male_homvars = hl.agg.count_where(mt.GT.is_hom_var() & (mt.is_female == False)),
        male_calls = hl.agg.count_where(hl.is_defined(mt.GT) & (mt.is_female == False)),
        female_hets = hl.agg.count_where(mt.GT.is_het() & (mt.is_female == True)),
        female_homvars = hl.agg.count_where(mt.GT.is_hom_var() & (mt.is_female == True)),
        female_calls = hl.agg.count_where(hl.is_defined(mt.GT) & (mt.is_female == True))
    )
    
    mt = mt.annotate_rows(
        call_rate = ( hl.case()
            .when(mt.locus.in_y_nonpar(), (mt.male_calls / num_males))
            .when(mt.locus.in_x_nonpar(), (mt.male_calls + 2*mt.female_calls) / (num_males + 2*num_females))
            .default((mt.male_calls + mt.female_calls) / (num_males + num_females)) ),
        AC = ( hl.case()
            .when(mt.locus.in_y_nonpar(), mt.male_homvars)
            .when(mt.locus.in_x_nonpar(), mt.male_homvars + mt.female_hets + 2*mt.female_homvars)
            .default(mt.male_hets + 2*mt.male_homvars + mt.female_hets + 2*mt.female_homvars) ),
        AN = ( hl.case()
            .when(mt.locus.in_y_nonpar(), mt.male_calls)
            .when(mt.locus.in_x_nonpar(), mt.male_calls + 2*mt.female_calls)
            .default(2*mt.male_calls + 2*mt.female_calls) ),
        pHWE = ( hl.case() 
            .when(mt.locus.in_y_nonpar() | mt.locus.in_mito(), 1.0)
            .when(mt.locus.in_x_nonpar(), hl.hardy_weinberg_test(
                hl.int32(mt.female_calls - mt.female_hets - mt.female_homvars), 
                hl.int32(mt.female_hets), 
                hl.int32(mt.female_homvars)).p_value)
            .default(hl.hardy_weinberg_test(
                hl.int32(mt.male_calls+mt.female_calls - mt.male_hets-mt.female_hets - mt.male_homvars-mt.female_homvars), 
                hl.int32(mt.male_hets + mt.female_hets), 
                hl.int32(mt.male_homvars + mt.female_homvars)).p_value) )
    )
    return mt

# Add annotations including pHWE
mt = sex_aware_variant_annotations_with_pHWE(mt)

# Filter non-variant sites (ought to be redundant with the way I structured things above)
mt = mt.filter_rows(mt.AC > 0, keep = True)

# Let's impose a modest variant call rate and pHWE filter
mt = mt.filter_rows((mt.call_rate < 0.8) | (mt.pHWE < 0.000000000001), keep = False)

# Define sex-aware sample call rate calculation
def sex_aware_sample_annotations(mt):
    num_y_non_par_vars = mt.aggregate_rows(hl.agg.count_where(mt.locus.in_y_nonpar()))
    num_all_other_vars = mt.aggregate_rows(hl.agg.count_where(~mt.locus.in_y_nonpar()))
    
    mt = mt.annotate_cols(
        sample_call_rate = ( hl.case()
            .when(mt.is_female == True, 
                hl.agg.count_where(hl.is_defined(mt.GT) & ~mt.locus.in_y_nonpar()) / num_all_other_vars)
            .default(
                hl.agg.count_where(hl.is_defined(mt.GT)) / (num_y_non_par_vars + num_all_other_vars)) )
    )
                        
    return mt

# Calculate sample call rates
mt = sex_aware_sample_annotations(mt)

# get sample-level stats to plot
hl.sample_qc(mt).cols().flatten().export(f"{cohort_prefix}_wes_final_annot_post_filter_qc_info.txt")


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

pedigree = hl.Pedigree.read(ped_uri)

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

de_novo_results.write(f"{cohort_prefix}_wes_final_denovo.ht", overwrite = True)

mt = mt.semi_join_rows(de_novo_results.key_by('locus', 'alleles'))

mt.rows().flatten().export(f"{cohort_prefix}_wes_final_denovo_vep.txt")

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
td.write(f"{cohort_prefix}_trio_tdt.mt", overwrite=True)

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
td.write(f"{cohort_prefix}_parent_aware_trio_tdt.mt", overwrite = True)


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