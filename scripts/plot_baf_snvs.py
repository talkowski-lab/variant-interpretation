import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import sys

merged_baf_tsv = sys.argv[1]

def human_size(interval_size):
    size_dict = {10**3: 'KB', 10**6: 'MB', 10**9: 'GB'}
    for unit, unit_str in size_dict.items():
        if interval_size / unit < 1:
            break
        prev_unit = unit
        prev_unit_str = unit_str
    new_interval_size = round(interval_size / prev_unit, 2)
    return f"{new_interval_size} {prev_unit_str}"

merged_baf = pd.read_csv(merged_baf_tsv, sep='\t')

try:
    merged_baf = merged_baf.reset_index()
    merged_baf.columns = ['locus', 'alleles', 'AB', 'locus_interval', 'SAMPLE', 'SV_type']
except:
    pass

merged_baf[['CHROM', 'POS']] = merged_baf.locus.str.split(':', expand=True)
merged_baf['chrom_int'] = merged_baf.CHROM.str.split('chr').str[1].replace({'X':23, 'Y':24}).astype(int)
merged_baf['POS'] = merged_baf.POS.astype(int)
merged_baf = merged_baf.sort_values(['chrom_int', 'POS'])

for locus_interval in merged_baf.locus_interval.unique():
    fig, ax = plt.subplots(figsize=(8, 4));
    locus_int_df = merged_baf[merged_baf.locus_interval==locus_interval]
    sample = locus_int_df.SAMPLE.unique()[0]
    start, end = (int(x) for x in locus_interval.split(':')[1].split('-'))
    interval_size = human_size(end - start)
    sv_type = locus_int_df.SV_type.unique()[0]
    median_baf = locus_int_df[(locus_int_df.AB<1)&(locus_int_df.AB>0)].AB.median()
    sns.scatterplot(data=locus_int_df, y='AB', x='POS', ax=ax);
    ax.set(ylim=(-0.1, 1.1), title=f"{locus_interval} ({interval_size} {sv_type}) in {sample}");
    ax.axhline(median_baf, linestyle='--', color='r');

    trans = transforms.blended_transform_factory(
        ax.get_yticklabels()[0].get_transform(), ax.transData)
    ax.text(1.03, median_baf, round(median_baf, 3), color="red", transform=trans, 
        ha="left", va="center")
    locus_str = locus_interval.replace(':', '_').replace('-', '_')
    plt.savefig(f"{locus_str}_{sv_type}_{sample}_baf.png");
