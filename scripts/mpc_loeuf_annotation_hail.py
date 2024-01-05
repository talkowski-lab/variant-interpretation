import hail as hl
import pandas as pd
import numpy as np
import sys
import os

vcf_metrics_tsv = sys.argv[1]
mpc_file = sys.argv[2]
loeuf_file = sys.argv[3]

mt = hl.import_table(vcf_metrics_tsv)

keys = hl.parse_variant(mt.ID, reference_genome='GRCh38')
mt = mt.annotate(alleles=keys.alleles, locus=keys.locus)

mpc = hl.read_table(mpc_file)
mpc = mpc.annotate(allele_array = [mpc.alleles[0], mpc.alleles[1]])
mpc = mpc.annotate(locus = mpc.locus_38)
mpc = mpc.key_by(mpc.locus, mpc.allele_array)

mpc_out = mpc[mt.locus, mt.alleles].MPC
mt = mt.annotate(MPC=mpc_out)

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