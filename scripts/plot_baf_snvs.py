import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import sys
import ast

merged_baf_tsv = sys.argv[1]
het_only = ast.literal_eval(sys.argv[2].capitalize())

def human_size(interval_size):
    size_dict = {10**3: 'KB', 10**6: 'MB', 10**9: 'GB'}
    for unit, unit_str in size_dict.items():
        if interval_size / unit < 1:
            break
        prev_unit = unit
        prev_unit_str = unit_str
    new_interval_size = round(interval_size / prev_unit, 2)
    return f"{new_interval_size} {prev_unit_str}"

def get_parental_origin(row):
    origin = 'unresolved'
    gt_sample, gt_mother, gt_father = row.GT_sample, row.GT_mother, row.GT_father
    if gt_sample=='0/1':
        if (gt_mother=='1/1') and (gt_father=='0/0'):
            origin = 'mother'
        elif (gt_father=='1/1') and (gt_mother=='0/0'):
            origin = 'father'
    return origin

merged_baf = pd.read_csv(merged_baf_tsv, sep='\t')

merged_baf[['CHROM', 'POS']] = merged_baf.locus.str.split(':', expand=True)
merged_baf['chrom_int'] = merged_baf.CHROM.str.split('chr').str[1].replace({'X':23, 'Y':24}).astype(int)
merged_baf['POS'] = merged_baf.POS.astype(int)
if het_only:
    merged_baf = merged_baf[merged_baf.GT_sample=='0/1']
merged_baf = merged_baf.sort_values(['chrom_int', 'POS'])
merged_baf['parental_origin'] = merged_baf.apply(get_parental_origin, axis=1).astype('category').cat.set_categories(['mother','father','unresolved'])    

for locus_interval in merged_baf.locus_interval.unique():
    fig, ax = plt.subplots(1, 3, figsize=(20, 5));
    locus_int_df = merged_baf[merged_baf.locus_interval==locus_interval]
    window_locus_interval = locus_int_df.window_locus_interval.unique()[0]
    sample = locus_int_df.SAMPLE.unique()[0]
    pipeline_id = locus_int_df.pipeline_id.unique()[0]
    start, end = (int(x) for x in locus_interval.split(':')[1].split('-'))
    interval_size = human_size(end - start)
    window_start, window_end = (int(x) for x in window_locus_interval.split(':')[1].split('-'))
    sv_type = locus_int_df.SV_type.unique()[0]

    fig.suptitle(f"{locus_interval} in {pipeline_id} ({interval_size} {sv_type})");
    for i, role in enumerate(['sample','father','mother']):
        median_baf_above = locus_int_df[(locus_int_df[f"AB_{role}"]<1)&(locus_int_df[f"AB_{role}"]>0.5)][f"AB_{role}"].median()
        median_baf_below = locus_int_df[(locus_int_df[f"AB_{role}"]>0)&(locus_int_df[f"AB_{role}"]<0.5)][f"AB_{role}"].median()
        sns.scatterplot(data=locus_int_df, y=f"AB_{role}", x='POS', ax=ax[i],  
                        hue='parental_origin' if i==0 else None, hue_order=['mother','father','unresolved'],
                        size='parental_origin' if i==0 else None, sizes={'unresolved': 9, 'mother': 12, 'father': 12} if i==0 else None,
                        s=9 if i!=0 else None,
                        palette={'mother':'mediumpurple', 'father':'mediumturquoise', 'unresolved':'silver'} if i==0 else None,
                        color='silver' if i!=0 else None);
        if i==0:
            ax[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1), ncol=3, fancybox=True, fontsize=8);
        ax[i].set(ylim=(-0.1, 1.1), title=f"AB_{role}");
        ax[i].axhline(median_baf_above, linestyle='--', color='indianred');
        ax[i].axhline(median_baf_below, linestyle='--', color='indianred');

        ax[i].axvline(start, linestyle='--', color='maroon');
        ax[i].axvline(end, linestyle='--', color='maroon');
        ax[i].axvspan(window_start, start, facecolor='0.2', alpha=0.02)
        ax[i].axvspan(end, window_end, facecolor='0.2', alpha=0.02)
        ax[i].margins(x=0)

        trans = transforms.blended_transform_factory(
            ax[i].get_yticklabels()[0].get_transform(), ax[i].transData)
        ax[i].text(1.03, median_baf_above, round(median_baf_above, 3), color="indianred", transform=trans, 
            ha="left", va="center");
        ax[i].text(1.03, median_baf_below, round(median_baf_below, 3), color="indianred", transform=trans, 
            ha="left", va="center");
    locus_str = locus_interval.replace(':', '_').replace('-', '_')
    plt.tight_layout();
    plt.savefig(f"{locus_str}_{sv_type}_{pipeline_id}.baf.png");
    plt.close();
