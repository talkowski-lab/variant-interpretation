from pyspark.sql import SparkSession
import hail as hl
import numpy as np
import sys

vcf_file = sys.argv[1]
vep_annotated_vcf_name = sys.argv[2]
cores = sys.argv[3]  # string
mem = int(np.floor(float(sys.argv[4])))

hl.init(spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g"})

header = hl.get_vcf_metadata(vcf_file) 
mt = hl.import_vcf(vcf_file, force_bgz=True, reference_genome='GRCh38')
mt = hl.vep(mt, config='vep_config.json', csq=True)
header['info']['CSQ'] = {'Description': hl.eval(mt.vep_csq_header), 'Number': '.', 'Type': 'String'}
mt = mt.annotate_rows(info = mt.info.annotate(CSQ=mt.vep))
hl.export_vcf(dataset=mt, output=vep_annotated_vcf_name, metadata=header)