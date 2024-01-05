import hail as hl
import pandas as pd
import numpy as np
import sys
import os

vcf_metrics_tsv = sys.argv[1]
mpc_dir = f"gs://{sys.argv[2]}"
mpc_chr22_file = sys.argv[3]
loeuf_file = sys.argv[4]

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
df.index = df.VarKey

df['gene_name'] = df.CSQ.str.split(',').str[0].str.split('|').str[3].replace({'': np.nan})
loeuf_vals = loeuf.loc[np.intersect1d(loeuf.index,df.gene_name.dropna()), 'LOEUF'].to_dict()
loeuf_tile_vals = loeuf.loc[np.intersect1d(loeuf.index,df.gene_name.dropna()), 'LOEUF_tile'].to_dict()
df['LOEUF'] = df.gene_name.map(loeuf_vals)
df['LOEUF_tile'] = df.gene_name.map(loeuf_tile_vals)

new_filename = os.path.basename(vcf_metrics_tsv).split('.tsv')[0] + '_with_mpc_loeuf.tsv'
df.to_csv(new_filename, sep='\t')