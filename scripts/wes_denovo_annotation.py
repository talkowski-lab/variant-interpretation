import hail as hl
import numpy as np
import sys

vcf_file = sys.argv[1]
cohort_prefix = sys.argv[2]
ped_uri = sys.argv[3]
gnomad_ht_uri = sys.argv[4]
mpc_dir = sys.argv[5]
mpc_chr22_file = sys.argv[6]
purcell5k = sys.argv[7]
cores = sys.argv[8]
mem = int(np.floor(float(sys.argv[9])))

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

mt = hl.import_vcf(vcf_file, reference_genome = 'GRCh38', array_elements_required=False, force_bgz=True)

# Step 1: Annotations

## Filter out discrepant samples from Somalier
ped = hl.import_table(ped_uri, impute=True).key_by('sample_id')
mt = mt.filter_cols(hl.is_defined(ped[mt.s]))

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

## Export mt
mt.write(f"{cohort_prefix}_wes_denovo_annot.mt", overwrite=True)