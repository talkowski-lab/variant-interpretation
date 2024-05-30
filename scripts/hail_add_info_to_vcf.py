import pandas as pd
import numpy as np
import hail as hl
import os
import sys

import gnomad
import gnomad.utils.vep

vcf_uri = sys.argv[1]
info_ht_uri = sys.argv[2]
vep_ht_uri = sys.argv[3]
qc_ht_uri = sys.argv[4]
cores = sys.argv[5]  # string
mem = int(np.floor(float(sys.argv[6])))

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")


mt = hl.import_vcf(vcf_uri, reference_genome='GRCh38', force_bgz=True, call_fields=[], array_elements_required=False)
mt = hl.split_multi_hts(mt)

header = hl.get_vcf_metadata(vcf_uri)
header['info']['CSQ'] = {'Description': gnomad.utils.vep.VEP_CSQ_HEADER, 'Number': '.', 'Type': 'String'}

# get row/INFO fields from INFO HT
info_ht = hl.read_table(info_ht_uri)
info_ht = info_ht.annotate( info = info_ht.info.annotate(
            QUALapprox = hl.case()
                .when(info_ht.info.QUALapprox > (2**31 - 1), (2**31 - 1))
                .default(info_ht.info.QUALapprox),
            AS_QUALapprox = hl.case()
                .when(info_ht.info.AS_QUALapprox > (2**31 - 1), (2**31 - 1))
                .default(info_ht.info.AS_QUALapprox)))
mt = mt.annotate_rows(info=info_ht[mt.row_key].info)
mt = mt.annotate_rows(qual=hl.float(info_ht[mt.row_key].info.AS_QUALapprox))

# remove all AC=0
mt = hl.variant_qc(mt)
mt = mt.filter_rows(mt.variant_qc.AC[1] > 0, keep = True)
mt = mt.annotate_rows(info=mt.info.annotate(AC=mt.variant_qc.AC[1:], 
                                    AF=mt.variant_qc.AF[1:],
                                    AN=mt.variant_qc.AN))
mt = mt.drop('variant_qc')

# get QUAL/FILTER info from QC HT
qc_ht = hl.read_table(qc_ht_uri)
mt = mt.annotate_rows(filters=qc_ht[mt.row_key].filters)

# get VEP info
vep_ht = hl.read_table(vep_ht_uri)
mt = mt.annotate_rows(info=mt.info.annotate(CSQ=gnomad.utils.vep.vep_struct_to_csq(vep_ht[mt.row_key].vep)))

output_filename = os.path.basename(vcf_uri).split('.vcf.bgz')[0] + '_info.vcf.bgz'
hl.export_vcf(mt, output_filename, metadata=header, tabix=True)