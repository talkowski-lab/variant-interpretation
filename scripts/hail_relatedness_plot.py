import pandas as pd
import numpy as np
import hail as hl
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
import ast
import datetime

rel_tsv = sys.argv[1]
cohort_prefix = sys.argv[2]
ped_uri = sys.argv[3]
cores = sys.argv[4]  # string
mem = int(np.floor(float(sys.argv[5])))

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

rel_df = pd.read_csv(rel_tsv, sep='\t')

fig, ax = plt.subplots(1, 2, figsize=(12, 5));
fig.suptitle(cohort_prefix);

sns.scatterplot(rel_df, x='ibd0', y='kin', hue='relationship', s=16, ax=ax[0],
               hue_order=['parent-child', 'siblings', 'second degree relatives', 'duplicate/twins', 'ambiguous', 'unrelated'], 
               palette={'parent-child': 'mediumpurple', 'siblings': 'mediumseagreen', 'second degree relatives': 'skyblue', 
                        'duplicate/twins': 'indianred', 'ambiguous': 'sandybrown', 'unrelated': 'silver'});
ax[0].set_title("Inferred relationship from VCF");

sns.scatterplot(rel_df, x='ibd0', y='kin', hue='ped_relationship', s=16, ax=ax[1],
               hue_order=['parent-child', 'siblings', 'related_other', 'unrelated'],
               palette={'parent-child': 'mediumpurple', 'siblings': 'mediumseagreen', 'related_other': 'skyblue', 'unrelated': 'silver'});
ax[1].set_title("Relationship from pedigree");

plt.tight_layout();
plt.savefig(f"{cohort_prefix}_relatedness_ibd0_kinship.png");