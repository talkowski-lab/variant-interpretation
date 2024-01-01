from pyspark.sql import SparkSession
import hail as hl
import numpy as np
import sys

vcf_file = sys.argv[1]
vep_annotated_vcf_name = sys.argv[2]
cores = sys.argv[3]  # string
mem = int(np.floor(float(sys.argv[4])))

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.memory": f"{mem}g"}, tmp_dir="tmp", local_tmpdir="tmp")

#split-multi
def split_multi_ssc(mt):
    mt = mt.annotate_rows(num_alleles = mt.alleles.size() ) # Add number of alleles at site before split
    # Now split
    sm = hl.split_multi(mt)
    pl = hl.or_missing(hl.is_defined(sm.PL),
                      (hl.range(0, 3).map(lambda i: hl.min(hl.range(0, hl.len(sm.PL))
       .filter(lambda j: hl.downcode(hl.unphased_diploid_gt_index_call(j), sm.a_index) == hl.unphased_diploid_gt_index_call(i))
       .map(lambda j: sm.PL[j])))))
    split_ds = sm.annotate_entries(GT = hl.downcode(sm.GT, sm.a_index),
                                   AD = hl.or_missing(hl.is_defined(sm.AD), [hl.sum(sm.AD) - sm.AD[sm.a_index], sm.AD[sm.a_index]]),
                                   PL = pl) 
        #GQ = hl.cond(hl.is_defined(pl[0]) & hl.is_defined(pl[1]) & hl.is_defined(pl[2]), hl.gq_from_pl(pl), sm.GQ) )
    mt = split_ds.drop('old_locus', 'old_alleles')
    return mt

header = hl.get_vcf_metadata(vcf_file) 
mt = hl.import_vcf(vcf_file, force_bgz=True, reference_genome='GRCh38')
mt = split_multi_ssc(mt)
# annotate cohort ac to INFO field (after splitting multiallelic)
mt = mt.annotate_rows(info=mt.info.annotate(cohort_AC=mt.info.AC[mt.a_index - 1],
                                           cohort_AF=mt.info.AF[mt.a_index - 1]))

mt = hl.vep(mt, config='vep_config.json', csq=True)
header['info']['CSQ'] = {'Description': hl.eval(mt.vep_csq_header), 'Number': '.', 'Type': 'String'}
mt = mt.annotate_rows(info = mt.info.annotate(CSQ=mt.vep))
hl.export_vcf(dataset=mt, output=vep_annotated_vcf_name, metadata=header)