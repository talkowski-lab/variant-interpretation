import hail as hl
import pandas as pd
import numpy as np
import sys
import os

vcf_metrics_tsv = sys.argv[1]
mpc_dir = sys.argv[2]
mpc_chr22_file = sys.argv[3]
loeuf_file = sys.argv[4]
cores = sys.argv[5]  # string
mem = int(np.floor(float(sys.argv[6])))

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

mt = hl.import_table(vcf_metrics_tsv)
keys = hl.parse_variant(mt.ID, reference_genome='GRCh38')
mt = mt.annotate(alleles=keys.alleles, locus=keys.locus)

# MPC
mpc = hl.read_table(mpc_dir)

mpc = mpc.annotate(allele_array = [mpc.alleles[0], mpc.alleles[1]])
mpc = mpc.annotate(locus = mpc.locus_38)
mpc = mpc.key_by(mpc.locus, mpc.allele_array)

mpc_chr22 = hl.import_table(mpc_chr22_file, types={"locus": hl.tlocus("GRCh38"), "alleles": hl.tarray(hl.tstr)})
mpc_chr22 = mpc_chr22.annotate(allele_array = [mpc_chr22.alleles[0], mpc_chr22.alleles[1]])
mpc_chr22 = mpc_chr22.key_by(mpc_chr22.locus, mpc_chr22.allele_array)
mpc_chr22 = mpc_chr22.annotate(MPC = hl.float64(mpc_chr22.MPC))

merged_mpc = mpc.select(mpc.MPC).union(mpc_chr22.select(mpc_chr22.MPC))

mt = mt.annotate(MPC=merged_mpc[mt.locus, mt.alleles].MPC)

# LOEUF
loeuf = pd.read_csv(loeuf_file, sep='\t')
loeuf.index = loeuf.gene_name

df = mt.to_pandas()

def get_genes_csq(csq):
    genes = []
    for ind_csq in csq:
        if ind_csq=='.':
            continue
        gene = ind_csq.split('|')[3]
        if gene!='':
            genes.append(gene)
    return list(set(genes))

if 'genes' not in df.columns:
    df['genes'] = df.CSQ.str.split(',').apply(get_genes_csq)

all_genes = df.genes.apply(pd.Series).stack().unique()

loeuf_vals = loeuf.loc[np.intersect1d(loeuf.index, all_genes), 'LOEUF'].to_dict()
loeuf_tile_vals = loeuf.loc[np.intersect1d(loeuf.index, all_genes), 'LOEUF_tile'].to_dict()

df['LOEUF'] = df.genes.apply(lambda gene_list: pd.Series(gene_list).map(loeuf_vals).min())
df['LOEUF_tile'] = df.genes.apply(lambda gene_list: pd.Series(gene_list).map(loeuf_tile_vals).min())

new_filename = os.path.basename(vcf_metrics_tsv).split('.tsv')[0] + '_with_mpc_loeuf.tsv.gz'
df.to_csv(new_filename, sep='\t', index=False)