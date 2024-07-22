import datetime
import pandas as pd
import hail as hl
import numpy as np
import sys
import os
import ast
import onnx
from gnomad.sample_qc.ancestry import apply_onnx_classification_model, apply_sklearn_classification_model, assign_population_pcs
from gnomad.utils.filtering import filter_to_adj

import matplotlib.pyplot as plt
import seaborn as sns

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
use_gnomad_rf = ast.literal_eval(sys.argv[11].capitalize())

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

mt = hl.import_vcf(vcf_uri, reference_genome='GRCh38', force_bgz=True, call_fields=[], array_elements_required=False)    
gnomad_mt = hl.import_vcf(gnomad_vcf_uri, reference_genome='GRCh38', force_bgz=True, call_fields=[], array_elements_required=False)    

loading_ht = hl.read_table(gnomad_loading_ht)
with hl.hadoop_open(gnomad_rf_onnx, "rb") as f:
    onx_fit = onnx.load(f)

# Filter the MT to high quality genotypes (GQ >= 20; DP >= 10; and AB >= 0.2 for het calls)
mt = filter_to_adj(mt)
common_entry_fields = [x for x in list(np.intersect1d(list(gnomad_mt.entry), list(mt.entry))) if x!='PGT']
mt = mt.select_entries(*common_entry_fields).union_cols(gnomad_mt.select_entries(*common_entry_fields), row_join_type='inner')

pop_labels_ht = hl.import_table(pop_labels_tsv)
pop_labels_ht = pop_labels_ht.annotate(s=pop_labels_ht.Sample).key_by('s')

pop_labels_ht = pop_labels_ht.annotate(SuperPop=pop_labels_ht.SuperPop.replace('CSA', 'SAS').replace('EUR', 'NFE'))
pop_labels_ht = pop_labels_ht.filter(pop_labels_ht.SuperPop!='OCE')

gnomad_pcs_ht = hl.experimental.pc_project(
    mt.GT, loading_ht.loadings, loading_ht.pca_af,
)

gnomad_pcs_ht = gnomad_pcs_ht.annotate(known_pop=pop_labels_ht[gnomad_pcs_ht.key].SuperPop.lower())

ht, model = assign_population_pcs(
    gnomad_pcs_ht,
    pc_cols=gnomad_pcs_ht.scores[:num_pcs],
    known_col='known_pop',
    fit=onx_fit if use_gnomad_rf else None,
    min_prob=min_prob,
    apply_model_func = apply_onnx_classification_model if use_gnomad_rf else apply_sklearn_classification_model
)

ht = ht.annotate(known_pop=gnomad_pcs_ht[ht.key].known_pop)
ancestry_df = ht.to_pandas()
ancestry_df.to_csv(cohort_prefix + '_inferred_ancestry.tsv', sep='\t', index=False)

# PLOT
ancestry_df.index = ancestry_df.s
pc_df = ancestry_df.pca_scores.apply(pd.Series)

ancestry_df = pd.concat([ancestry_df, pc_df.rename({i: f"PC_{i}" for i in pc_df.columns}, axis=1)], axis=1)

POP_COLORS = {
    "afr": "#941494",
    "ami": "#FFC0CB",
    "amr": "#ED1E24",
    "asj": "#FF7F50",
    "eas": "#108C44",
    "eur": "#6AA5CD",
    "fin": "#002F6C",
    "mid": "#33CC33",
    "nfe": "#6AA5CD",
    "oth": "#ABB9B9",
    "unknown": "#ABB9B9",
    "sas": "#FF9912"}

fig, ax = plt.subplots(figsize=(6, 5));
sns.scatterplot(data=ancestry_df, x='PC_1', y='PC_2', hue='known_pop', palette=POP_COLORS, ax=ax, alpha=0.1, s=70);
sns.scatterplot(data=ancestry_df[ancestry_df.known_pop.isna()], x='PC_1', y='PC_2', hue='pop', palette=POP_COLORS, ax=ax, legend=False);
leg = plt.legend();
for lh in leg.legend_handles:
    lh.set_alpha(1);
ax.set(title=cohort_prefix);
plt.tight_layout();
plt.savefig(f"{cohort_prefix}_inferred_ancestry.png");