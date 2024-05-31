import datetime
import pandas as pd
import hail as hl
import numpy as np
import sys
import ast
import os

file = sys.argv[1]
cohort_prefix = sys.argv[2]
ped_uri = sys.argv[3]
gnomad_ht_uri = sys.argv[4]
mpc_ht_uri = sys.argv[5]
cores = sys.argv[6]
mem = int(np.floor(float(sys.argv[7])))
bucket_id = sys.argv[8]
genome_build = sys.argv[9]

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

if file.split('.')[-1] == 'mt':
    mt = hl.read_matrix_table(file)
    prefix = os.path.basename(file).split('.mt')[0]
else:
    mt = hl.import_vcf(file, reference_genome = genome_build, array_elements_required=False, force_bgz=True)
    prefix = os.path.basename(file).split('.vcf')[0]

# Step 1: Annotations

## gnomAD exome frequency annotations
gnomad_ht = hl.read_table(gnomad_ht_uri)

# annotate with gnomAD exome frequencies; use "non-neuro" allele frequencies 
mt = mt.annotate_rows(gnomad_non_neuro_AF = 
                      gnomad_ht.index(mt.row_key).freq[hl.eval(gnomad_ht.freq_index_dict["non_neuro"])].AF)

## MPC annotations
mpc = hl.read_table(mpc_ht_uri).key_by('locus','alleles')
mt = mt.annotate_rows(MPC=mpc[mt.locus, mt.alleles].mpc)

## pAB annotations

# Add pAB entry field for downstream filtering -- note, only for hets
mt = mt.annotate_entries(AB = mt.AD[1] / hl.sum(mt.AD), 
                         pAB = hl.or_missing(mt.GT.is_het(), 
                                             hl.binom_test(mt.AD[1], hl.sum(mt.AD), 0.5, 'two-sided')))

# Run and export Hail's built-in sample QC function
qc_filename = f"{prefix}_wes_post_annot_sample_QC_info.txt"
hl.sample_qc(mt, 'sample_qc').cols().flatten().export(
    qc_filename)

pd.Series([qc_filename]).to_csv('qc_out.txt',index=False, header=None)

## Export mt
filename = f"{bucket_id}/hail/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))}/{prefix}_wes_denovo_annot.mt"
pd.Series([filename]).to_csv('mt_uri.txt',index=False, header=None)

mt.write(filename, overwrite=True)
