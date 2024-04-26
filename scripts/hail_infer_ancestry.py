import datetime
import pandas as pd
import hail as hl
import numpy as np
import sys
import os
import onnx
from gnomad.sample_qc.ancestry import apply_onnx_classification_model, assign_population_pcs
from gnomad.utils.filtering import filter_to_adj

vcf_uri = sys.argv[1]
gnomad_vcf_uri = sys.argv[2]
gnomad_loading_ht = sys.argv[3]
gnomad_rf_onnx = sys.argv[4]
pop_labels_tsv = sys.argv[5]
num_pcs = int(sys.argv[6])
min_prob = float(sys.argv[7])
cohort_prefix = sys.argv[8]
cores = sys.argv[9]
mem = int(np.floor(float(sys.argv[10])))

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

mt = hl.import_vcf(vcf_uri, reference_genome='GRCh38', force_bgz=True, call_fields=[], array_elements_required=False)    
gnomad_mt = hl.import_vcf(gnomad_vcf_uri, reference_genome='GRCh38', force_bgz=True, call_fields=[], array_elements_required=False)    
loading_ht = hl.read_table(gnomad_loading_ht)
with hl.hadoop_open(gnomad_rf_onnx, "rb") as f:
    onx_fit = onnx.load(f)

# Filter the MT to high quality genotypes (GQ >= 20; DP >= 10; and AB >= 0.2 for het calls)
mt = filter_to_adj(mt)
common_entry_fields = [x for x in list(np.intersect1d(list(gnomad_mt.entry), list(mt.entry))) if x!='PGT']
mt = mt.select_entries(*common_entry_fields).union_cols(gnomad_mt.select_entries(*common_entry_fields), row_join_type='outer')

pop_labels_ht = hl.import_table(pop_labels_tsv)
pop_labels_ht = pop_labels_ht.annotate(s=pop_labels_ht.Sample).key_by('s')

gnomad_pcs_ht = hl.experimental.pc_project(
    mt.GT, loading_ht.loadings, loading_ht.pca_af,
)

gnomad_pcs_ht = gnomad_pcs_ht.annotate(known_pop=pop_labels_ht[gnomad_pcs_ht.key].SuperPop.lower())

ht, model = assign_population_pcs(
    gnomad_pcs_ht,
    pc_cols=gnomad_pcs_ht.scores[:num_pcs],
    known_col='known_pop',
    fit=onx_fit,
    min_prob=min_prob,
    apply_model_func = apply_onnx_classification_model,
)

ht = ht.annotate(known_pop=gnomad_pcs_ht[ht.key].known_pop)
df = ht.to_pandas()
df.to_csv(cohort_prefix + '_inferred_ancestry.tsv', sep='\t', index=False)