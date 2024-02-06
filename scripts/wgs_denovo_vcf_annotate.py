import hail as hl
import pandas as pd
import numpy as np
import sys
import ast
import os

vcf_file = sys.argv[1]
vep_uri = sys.argv[2]
mpc_dir = sys.argv[3]
mpc_chr22_file = sys.argv[4]
loeuf_file = sys.argv[5]
header_file = sys.argv[6]
file_ext = sys.argv[7]
sample = sys.argv[8]
cores = sys.argv[9]
mem = int(np.floor(float(sys.argv[10])))

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

mt = hl.import_vcf(vcf_file, force_bgz=True, array_elements_required=False, call_fields=[], header_file=header_file, reference_genome='GRCh38')

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

mt = mt.annotate_rows(info = mt.info.annotate(MPC=merged_mpc[mt.locus, mt.alleles].MPC))

mt.info.CSQ.export(f"{os.path.basename(vcf_file).split(file_ext)[0]}_CSQ.txt")

mt.info.MPC.export(f"{os.path.basename(vcf_file).split(file_ext)[0]}_MPC.txt")

df = pd.concat([pd.read_csv(f"{os.path.basename(vcf_file).split(file_ext)[0]}_CSQ.txt", sep='\t'),
                pd.read_csv(f"{os.path.basename(vcf_file).split(file_ext)[0]}_MPC.txt", sep='\t').iloc[:,2]], axis=1)

df.columns = ['locus', 'alleles', 'CSQ', 'MPC']
        
metadata = hl.get_vcf_metadata(vep_uri)
csq_columns = metadata['info']['CSQ']['Description'].split('Format: ')[1].split('|')

loeuf = pd.read_csv(loeuf_file, sep='\t')
loeuf.index = loeuf.gene_name

df['CSQ'] = df.CSQ.replace({np.nan, '[]'}).apply(ast.literal_eval)

def combine_csq(csq, col_num):
    csqs = []
    for ind_csq in csq:
        csqs.append(ind_csq.split('|')[col_num])
    return csqs

df['TYPE'] = df.CSQ.apply(combine_csq, args=[csq_columns.index('VARIANT_CLASS')]).apply(lambda lst: pd.Series(lst).unique()[0]).map({'insertion': 'Indel', 'deletion': 'Indel'})
df['all_genes'] = df.CSQ.apply(combine_csq, args=[csq_columns.index('SYMBOL')]).apply(lambda lst: pd.Series(lst).unique())

all_genes = df.all_genes.apply(pd.Series).stack().unique()

loeuf_vals = loeuf.loc[np.intersect1d(loeuf.index, all_genes), 'LOEUF'].to_dict()
loeuf_tile_vals = loeuf.loc[np.intersect1d(loeuf.index, all_genes), 'LOEUF_tile'].to_dict()

df['LOEUF'] = df.all_genes.apply(lambda gene_list: pd.Series(gene_list).map(loeuf_vals).min())
df['LOEUF_tile'] = df.all_genes.apply(lambda gene_list: pd.Series(gene_list).map(loeuf_tile_vals).min())

if sample!='false':
    df['SAMPLE'] = sample

df.to_csv(f"{os.path.basename(vcf_file).split(file_ext)[0]}_annot.tsv", sep='\t', index=False)