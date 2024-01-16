import sys
import hail as hl
import numpy as np 

vcf_file = sys.argv[1]
records_per_shard = sys.argv[2]
prefix = sys.argv[3]
cores = sys.argv[4]
mem = int(np.floor(float(sys.argv[5])))

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

mt = hl.import_vcf(vcf_file, force_bgz=True, reference_genome='GRCh38')

# num_variants = mt.count_rows()

i=0
# while True:
#     if i*records_per_shard >= num_variants:
#         break
#     if ((i+1)%50==0):
#         print(f"Doing shard {i+1}...")
#     mt_shard = mt[i*records_per_shard:(i+1)*records_per_shard, :]
#     hl.export_vcf(mt_shard, f"{prefix}.shard_{i}.vcf.gz")
#     i+=1

while True:
    if ((i+1)%50==0):
        print(f"Doing shard {i+1}...")
    try:
        mt_shard = mt.head(records_per_shard)
        mt = mt.tail()
        hl.export_vcf(mt_shard, f"{prefix}.shard_{i}.vcf.gz")
        i+=1
    except:
        break

print(f"VCF split into {i+1} shards.")
