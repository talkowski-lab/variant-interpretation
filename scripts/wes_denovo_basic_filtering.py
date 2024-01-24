import hail as hl
import numpy as np
import sys

annot_mt = sys.argv[1]
cohort_prefix = sys.argv[2]
corrected_ped = sys.argv[3]
bucket_id = sys.argv[4]
cores = sys.argv[5]
mem = int(np.floor(float(sys.argv[6])))

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

# Step 2: Filtering

## Basic genotype filtering

# Following filters:

# - All calls: Depth below 7 or above 1000
# - HomRef calls: GQ below 25
# - HomVar calls: PL(HomRef) below 25
# - Het calls:PL(HomRef) below 25

mt = hl.read_matrix_table(annot_mt)

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

ped = hl.import_table(corrected_ped, impute=True)

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

# export mt
mt.write(f"{bucket_id}/hail/{cohort_prefix}_wes_denovo_basic_filtering.mt")