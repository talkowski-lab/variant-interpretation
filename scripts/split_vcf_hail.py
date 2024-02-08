import sys
import hail as hl
import numpy as np 

file = sys.argv[1]
n_shards = int(sys.argv[2])
prefix = sys.argv[3]
cores = sys.argv[4]
mem = int(np.floor(float(sys.argv[5])))

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

is_vcf = file.split('.')[-1] != 'mt'
if not is_vcf:
    mt = hl.read_matrix_table(file)
else:
    mt = hl.import_vcf(file, force_bgz=True, array_elements_required=False, reference_genome='GRCh38')
    header = hl.get_vcf_metadata(file) 

# for haploid (e.g. chrY)
mt = mt.annotate_entries(
    GT = hl.if_else(
             mt.GT.ploidy == 1, 
             hl.call(mt.GT[0], mt.GT[0]),
             mt.GT)
)

if n_shards!=0:
    mt = mt.repartition(n_shards)

hl.export_vcf(mt, output=prefix+'.vcf.bgz', parallel='header_per_shard', metadata=header if is_vcf else None)