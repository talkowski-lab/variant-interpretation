import sys
import hail as hl
import numpy as np 

vcf_file = sys.argv[1]
n_shards = int(sys.argv[2])
prefix = sys.argv[3]
cores = sys.argv[4]
mem = int(np.floor(float(sys.argv[5])))

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

mt = hl.import_vcf(vcf_file, force_bgz=True, reference_genome='GRCh38')
header = hl.get_vcf_metadata(vcf_file) 

mt = mt.repartition(n_shards)
hl.export_vcf(mt, output=prefix+'.vcf.bgz', parallel='header_per_shard', metadata=header)