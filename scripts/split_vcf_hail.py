import sys
import hail as hl
import numpy as np 

file = sys.argv[1]
n_shards = int(sys.argv[2])
records_per_shard = int(sys.argv[3])
prefix = sys.argv[4]
cores = sys.argv[5]
mem = int(np.floor(float(sys.argv[6])))
build = sys.argv[7]
# row_fields_to_keep = sys.argv[7].split(',')

# if row_fields_to_keep[0] == 'false':
#     row_fields_to_keep = []

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

is_vcf = file.split('.')[-1] != 'mt'
if not is_vcf:
    mt = hl.read_matrix_table(file, n_partitions=n_shards if n_shards!=0 else None)
else:
    mt = hl.import_vcf(file, force_bgz=True, array_elements_required=False, call_fields=[], reference_genome=build, find_replace=('nan', '.'),
                       n_partitions=n_shards if n_shards!=0 else None)
    header = hl.get_vcf_metadata(file) 

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
    mt = mt.repartition(n_shards)

# put all in INFO to be kept when exported to VCF--this doesn't really work because of header/metadata issues...
# for field in row_fields_to_keep:
#     mt = mt = mt.annotate_rows(info = mt.info.annotate(**{field: getattr(mt, field)}))
    
tot_num_records = mt.count_rows()
if tot_num_records > 0:
    hl.export_vcf(mt, output=prefix+'.vcf.bgz', parallel='header_per_shard', metadata=header if is_vcf else None)