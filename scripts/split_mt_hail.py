import sys
import pandas as pd
import hail as hl
import numpy as np 
import datetime
import os

mt_uri = sys.argv[1]
n_shards = int(sys.argv[2])
records_per_shard = int(sys.argv[3])
cores = sys.argv[4]
mem = int(np.floor(float(sys.argv[5])))
bucket_id = sys.argv[6]

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

mt = hl.read_matrix_table(mt_uri)

try:
    # for haploid (e.g. chrY)
    mt = mt.annotate_entries(
        GT = hl.if_else(
                mt.GT.ploidy == 1, 
                hl.call(mt.GT[0], mt.GT[0]),
                mt.GT)
    )
except:
    pass

if records_per_shard!=0:
    tot_num_records = mt.count_rows()
    n_shards = int(np.ceil(tot_num_records / records_per_shard))

if n_shards!=0:
    mt = mt.repartition(n_shards)

mt_dirs = []
for i in range(n_shards):
    print(f"writing shard {i}/{n_shards}...")
    dir_name = f"{bucket_id}/hail/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))}/{os.path.basename(mt_uri).split('.mt')[0]}_shard_{i}.mt"
    mt_dirs.append(dir_name)
    mt._filter_partitions([i]).write(dir_name, overwrite=True)

pd.Series([mt_dirs]).to_csv('mt_uris.txt',index=False, header=None)
