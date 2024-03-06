

# TODO: test in notebook
import datetime
import pandas as pd
import hail as hl
import numpy as np
import sys

file = sys.argv[1]
cohort_prefix = sys.argv[2]
ped_uri = sys.argv[3]
gnomad_ht_uri = sys.argv[4]
mpc_dir = sys.argv[5]
mpc_chr22_file = sys.argv[6]
purcell5k = sys.argv[7]
cores = sys.argv[8]
mem = int(np.floor(float(sys.argv[9])))
bucket_id = sys.argv[10]

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

if file.split('.')[-1] == 'mt':
    mt = hl.read_matrix_table(file)
else:
    mt = hl.import_vcf(file, reference_genome = 'GRCh38', array_elements_required=False, force_bgz=True)


## PCA on 5k variants only
# Create filtered mt
# Load lifted-over Purcell 5k interval list
p5k = hl.import_locus_intervals(purcell5k, 
                                 reference_genome='GRCh38') #few variants that are likely most useful (PCA and relatedness)
mt5k = mt.filter_rows(hl.is_defined(p5k[mt.locus]), keep = True)

#5k variants only
eigenvalues_5k, score_table_5k, loading_table_5k = hl.hwe_normalized_pca(mt5k.GT, k=20, compute_loadings=True)

if bucket_id == 'false':
    score_table_5k.write(f"{cohort_prefix}_wes_pca_score_table_5k.ht", overwrite = True)
    loading_table_5k.write(f"{cohort_prefix}_wes_pca_loading_table_5k.ht", overwrite = True)

else:
    score_table_file = f"{bucket_id}/hail/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))}/{cohort_prefix}_wes_pca_score_table_5k.ht"
    loading_table_file = f"{bucket_id}/hail/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))}/{cohort_prefix}_wes_pca_loading_table_5k.ht"
        
    mt_uris = [score_table_file, loading_table_file]
    score_table_5k.write(score_table_file, overwrite=True)
    loading_table_5k.write(loading_table_file, overwrite=True)


