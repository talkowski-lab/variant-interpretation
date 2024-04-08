# from wgs_bagging_pu_cross_validation.ipynb
import gzip
import os
import glob
import io
import ast
import warnings
import seaborn as sns
import pandas as pd
import numpy as np
import xgboost as xgb
import matplotlib.pyplot as plt
import sklearn.metrics
import sklearn.model_selection

from baggingPU import BaggingClassifierPU
from sklearn.metrics import RocCurveDisplay
from google.cloud import storage
from sklearn.ensemble import RandomForestClassifier
from io import StringIO

plt.rcParams.update({'font.size': 10})

import sys
import pandas as pd
import numpy as np
import sklearn.metrics
import sklearn.model_selection
import ast
import pickle 

from sklearn.ensemble import RandomForestClassifier

vcf_metrics_tsv = sys.argv[1]
ultra_rare_variants_tsv = sys.argv[2]
cohort_prefix = sys.argv[3]
var_type = sys.argv[4]
variant_features = sys.argv[5].split(',')
sample_features = sys.argv[6].split(',')
vqslod_cutoff = float(sys.argv[7])
prop_dn = float(sys.argv[8])
ultra_rare_rep_regions = sys.argv[9]
rep_regions = sys.argv[10]
known_vars_uri = sys.argv[11]
metric = sys.argv[12]  # ['roc-auc', 'accuracy', 'f1', 'fp_fn_ratio']
n_estimators_rf = int(sys.argv[13])
n_bags = int(sys.argv[14])

def fp_fn_ratio(y, y_pred):
    FP = ((y==0) & (y_pred==1)).sum()
    TP = ((y==1) & (y_pred==1)).sum()
    FN = ((y==1) & (y_pred==0)).sum()
    TN = ((y==0) & (y_pred==0)).sum()
    FN_rate = FN / (FN + TN)
    FP_rate = FP / (FP + TP)
    if np.isnan(FP_rate):
        FP_rate = FN_rate
    return float(FP_rate / FN_rate)

metrics_to_funcs = {'roc-auc': sklearn.metrics.roc_auc_score, 'accuracy': True, 'f1': sklearn.metrics.f1_score,
                    'fp_fn_ratio': fp_fn_ratio}

known_vars_exist = (known_vars_uri!='false')

# variable params
# n_bags, numeric, p_thr, n_estimators_rf, metric
# choose best n_estimators_rf for metric
# choose n_bags --> plot (likely to not have much of an effect)

# https://heartbeat.comet.ml/positive-and-unlabelled-learning-recovering-labels-for-data-using-machine-learning-59c1def5452f
# https://roywrightme.wordpress.com/2017/11/16/positive-unlabeled-learning/

# Functions

def load_variants(vcf_metrics_tsv, ultra_rare_variants_tsv, var_type): 
    ultra_rare = pd.read_csv(ultra_rare_variants_tsv, sep='\t')
    ultra_rare['Indel_type'] = ultra_rare.apply(lambda x: 'Insertion' if (len(x.ALT) - len(x.REF)) > 0 else 'Deletion', axis=1)
    ultra_rare.loc[ultra_rare.TYPE=='SNV', 'Indel_type'] = 'SNV'
    ultra_rare = ultra_rare[ultra_rare.TYPE==var_type]
    # ultra_rare = ultra_rare[~ultra_rare.SAMPLE.isin(outlier_samples)]

    final_output = pd.read_csv(vcf_metrics_tsv, sep='\t')
    final_output['Indel_type'] = final_output.apply(lambda x: 'Insertion' if (len(x.ALT) - len(x.REF)) > 0 else 'Deletion', axis=1)
    final_output.loc[final_output.TYPE=='SNV', 'Indel_type'] = 'SNV'
    final_output = final_output[final_output.TYPE==var_type]
    # final_output = final_output[~final_output.SAMPLE.isin(outlier_samples)]

    ultra_rare = ultra_rare[(~ultra_rare.AD_father.isna()) & (~ultra_rare.AD_mother.isna())].reset_index(drop=True)

    ultra_rare['VarKey'] = ultra_rare.ID + ':' + ultra_rare.SAMPLE

    # ultra_rare['AB_mother'] = ultra_rare.AD_mother.apply(ast.literal_eval).apply(min) / ultra_rare.AD_mother.apply(ast.literal_eval).apply(sum)
    # ultra_rare['AB_father'] = ultra_rare.AD_father.apply(ast.literal_eval).apply(min) / ultra_rare.AD_father.apply(ast.literal_eval).apply(sum)
    # ultra_rare['AB_sample'] = ultra_rare.AD_sample.apply(ast.literal_eval).apply(min) / ultra_rare.AD_sample.apply(ast.literal_eval).apply(sum)

    ultra_rare = ultra_rare[~ultra_rare.PL_sample.isna()]
    ultra_rare[['PL_sample_0.0', 'PL_sample_0.1', 'PL_sample_1.1']] = ultra_rare.PL_sample.str.strip('][').str.split(', ', expand=True).astype(int)
    ultra_rare['GQ_mean'] = ultra_rare[['GQ_sample', 'GQ_mother', 'GQ_father']].mean(axis=1)

    # https://gatk.broadinstitute.org/hc/en-us/articles/360040507131-FilterVcf-Picard-
    # Allele balance is calculated for heterozygotes as the number of bases supporting the least-represented allele over the total number of base observations. Different heterozygote genotypes at the same locus are measured independently
    
    for samp in ['sample', 'mother', 'father']:
        final_output[f"DPC_{samp}"] = final_output[f"AD_{samp}"].apply(ast.literal_eval).apply(sum)
        final_output[f"AB_{samp}"] = final_output[f"AD_{samp}"].apply(ast.literal_eval).str[1] / final_output[f"DPC_{samp}"]

    # ultra_rare['GQ_parent'] = ultra_rare.apply(lambda row: row.GQ_mother if row.GT_mother in ['0/0', '0|0']
    #                                         else row.GQ_father, axis=1)
    # ultra_rare['AB_parent'] = ultra_rare.apply(lambda row: row.AB_mother if row.GT_mother in ['0/0', '0|0'] 
    #                                             else row.AB_father, axis=1)
    # ultra_rare['DPC_parent'] = ultra_rare.apply(lambda row: row.DPC_mother if row.GT_mother in ['0/0', '0|0'] 
    #                                             else row.DPC_father, axis=1)
    # ultra_rare['VAF_parent'] = ultra_rare.apply(lambda row: row.VAF_mother if row.GT_mother in ['0/0', '0|0'] 
    #                                             else row.VAF_father, axis=1)
    ultra_rare['GQ_parent'] = ultra_rare[['GQ_mother', 'GQ_father']].min(axis=1)
    ultra_rare['AB_parent'] = ultra_rare[['AB_mother', 'AB_father']].min(axis=1)
    ultra_rare['DPC_parent'] = ultra_rare[['DPC_mother', 'DPC_father']].min(axis=1)
    ultra_rare['VAF_parent'] = ultra_rare[['VAF_mother', 'VAF_father']].min(axis=1)
    
    final_output['GQ_parent'] = final_output[['GQ_mother', 'GQ_father']].min(axis=1)
    final_output['AB_parent'] = final_output[['AB_mother', 'AB_father']].min(axis=1)
    final_output['DPC_parent'] = final_output[['DPC_mother', 'DPC_father']].min(axis=1)
    final_output['VAF_parent'] = final_output[['VAF_mother', 'VAF_father']].min(axis=1)

    final_output_raw = final_output.copy();
    ultra_rare_raw = ultra_rare.copy();
    return final_output, ultra_rare, final_output_raw, ultra_rare_raw

def filter_variants(final_output, ultra_rare, final_output_raw, ultra_rare_raw):
    if (final_output['VQSLOD']!='.').sum()!=0:
        final_output = final_output[final_output['VQSLOD']!='.'].reset_index(drop=True)
        final_output['VQSLOD'] = final_output.VQSLOD.astype(float)

    final_output = final_output[final_output.VQSLOD>vqslod_cutoff]
    final_output = final_output[final_output.LEN<=50]
    final_output = final_output[final_output.FILTER=='PASS']

    # TODO: remove when filter-rare-variants-hail is fixed?
    # ultra_rare = ultra_rare[ultra_rare.GQ_sample>=99]
    # ultra_rare = ultra_rare[(ultra_rare.GQ_mother>=30)&(ultra_rare.GQ_father>=30)]

    ultra_rare = ultra_rare[ultra_rare.POLYX<=10]
    ultra_rare = ultra_rare[ultra_rare.LEN<=50]
    try:
        ultra_rare = ultra_rare[~ultra_rare.VQSLOD.isna()]
        ultra_rare = ultra_rare[ultra_rare.VQSLOD>vqslod_cutoff]
    except:
        pass
    ultra_rare = ultra_rare[(~ultra_rare.VAF_sample.isna())&
                            (~ultra_rare.VAF_mother.isna())&
                            (~ultra_rare.VAF_father.isna())].reset_index(drop=True)

    ultra_rare['label'] = 1  # positive
    final_output['label'] = 0
    # overlapping variants set to unlabeled
    ultra_rare = ultra_rare[~ultra_rare.VarKey.isin(np.intersect1d(ultra_rare.VarKey, final_output.VarKey))]
  
    merged_output = pd.concat([ultra_rare, final_output]).reset_index(drop=True)
    merged_output['var_type'] = merged_output.label.map({1: 'ultra-rare', 0: 'unlabeled'})

    return final_output, ultra_rare, merged_output

def BaggingPU(X, y, kf, model, n_bags, n_jobs=-1):
    y_pred_bag = np.zeros(y.size)
    output_bag = np.zeros(y.size)
    classifiers = []
    
    for i, (train_index, test_index) in enumerate(kf.split(X)):
        print(f"----------------- SPLIT {i+1}/{kf.get_n_splits()} -----------------")
        X_train = X.iloc[train_index,:]
        y_train = y.iloc[train_index]
        X_test = X.iloc[test_index,:]
        y_test = y.iloc[test_index]   
        bc = BaggingClassifierPU(
            model, 
            n_estimators = n_bags, n_jobs = n_jobs, 
            max_samples = sum(y_train==0)  # all unlabeled
        )
        bc.fit(X_train, y_train)
        y_pred_bag[test_index] = bc.predict(X_test)
        output_bag[test_index] = bc.predict_proba(X_test)[:,1]
        classifiers.append(bc)
    return y_pred_bag, output_bag, classifiers


def runBaggingPU_RF(X, y, model, merged_output, numeric, suffix, n_bags=10):
    kf = sklearn.model_selection.KFold(n_splits=5)
    y_pred_bag, output_bag, classifiers_bag = BaggingPU(X, y, kf, model, n_bags=n_bags, n_jobs=-1)

    print('---- {} ----'.format('Bagging PU'))
    print(sklearn.metrics.confusion_matrix(y, y_pred_bag))

    importances_bag = pd.DataFrame(np.array([[classifiers_bag[j].estimators_[i].feature_importances_ 
            for i in range(len(classifiers_bag[j].estimators_))] 
                for j in range(len(classifiers_bag))]).reshape((len(classifiers_bag)*len(classifiers_bag[0]),len(numeric))),
                    columns=numeric)

    results = pd.DataFrame({'label': y, 
                            'VarKey': merged_output.iloc[X.index].VarKey,
                            f'predict_proba_bag_{suffix}': output_bag,
                            f'pred_bag_{suffix}': y_pred_bag.astype(int)
                            })
    if known_vars_exist:
        results['is_known'] = merged_output.iloc[X.index].is_known.astype(str)
        results.loc[results.label==1, 'is_known'] = 'ultra-rare'    
    return results, classifiers_bag




def get_importances_oob_scores(X, y, merged_output, numeric, suffix, n_estimators_rf=100, n_bags=10, oob_score=True):
    model = RandomForestClassifier(warm_start=True, oob_score=oob_score, n_estimators=n_estimators_rf)
    results, estimators = runBaggingPU_RF(X, y, model, merged_output, numeric, suffix, n_bags=n_bags)

    importances = pd.DataFrame(np.array([[estimators[j].estimators_[i].feature_importances_ 
                for i in range(len(estimators[j].estimators_))] 
                    for j in range(len(estimators))]).reshape((len(estimators)*len(estimators[0]),len(numeric))),
                        columns=numeric)
    oob_scores = pd.DataFrame(np.array([[estimators[j].estimators_[i].oob_score_ 
                for i in range(len(estimators[j].estimators_))] 
                    for j in range(len(estimators))]).T)
    
    return results, estimators, importances, oob_scores


from sklearn.metrics import RocCurveDisplay

def run_baggingPU_level_features(X, y, merged_output, numeric,
                                 n_estimators_rf, n_bags, metric, suffix):
    
    results_optimized, estimators_optimized, importances_optimized, oob_scores_optimized = get_importances_oob_scores(X, y, merged_output, numeric, suffix,
                                                                                                                      n_estimators_rf=n_estimators_rf, n_bags=n_bags, oob_score=metrics_to_funcs[metric])

    opt_auc_scores = pd.DataFrame()

    if known_vars_exist:
#         fig, ax = plt.subplots(1, 2, figsize=(10,4));
#         fig.suptitle(f"{cohort_prefix} {var_type}");

        for p_thr in np.arange(0, 0.5, 0.05): 
            y_pred = [1 if x>p_thr else 0 for x in results_optimized[f'predict_proba_bag_{suffix}']]
#             RocCurveDisplay.from_predictions(y, y_pred, ax=ax[0], name=f"p_thr={round(p_thr, 2)}");
#             ax[0].set(title=f"unlabeled vs. ultra-rare");
            opt_auc_scores.loc[p_thr, 'roc_auc_score'] = sklearn.metrics.roc_auc_score(results_optimized.label, y_pred)

            y_pred_known = [1 if x>p_thr else 0 for x in results_optimized[results_optimized.is_known!='ultra-rare'][f'predict_proba_bag_{suffix}']]
            y_known = results_optimized[results_optimized.is_known!='ultra-rare'].is_known.apply(ast.literal_eval).astype(int)
#             RocCurveDisplay.from_predictions(y_known, y_pred_known, ax=ax[1], name=f"p_thr={round(p_thr, 2)}");
#             ax[1].set(title=f"known variants vs. unknown");

    else: 
#         fig, ax = plt.subplots(figsize=(4,4));
#         ax.set_title(f"{cohort_prefix} {var_type},\nunlabeled vs. ultra-rare");

        for p_thr in np.arange(0, 0.5, 0.05): 
            y_pred = [1 if x>p_thr else 0 for x in results_optimized[f'predict_proba_bag_{suffix}']]
#             RocCurveDisplay.from_predictions(y, y_pred, ax=ax, name=f"p_thr={round(p_thr, 2)}");
            opt_auc_scores.loc[p_thr, 'roc_auc_score'] = sklearn.metrics.roc_auc_score(results_optimized.label, y_pred)

#     plt.savefig(f"{cohort_prefix}_{var_type}_{metric}_RF_roc_auc_curve_p_thr.png");
#     plt.show();

    best_p_thr = round(opt_auc_scores.roc_auc_score.idxmax(), 2) 

    results_optimized[f'pred_bag_optimized_{suffix}'] = [1 if x>best_p_thr else 0 for x in results_optimized[f'predict_proba_bag_{suffix}']]
    results_optimized[f'features_{suffix}'] = ', '.join(numeric)
    results_optimized['metric'] = metric
    results_optimized[f'p_thr_{suffix}'] = best_p_thr

    return results_optimized, estimators_optimized, importances_optimized, oob_scores_optimized

# Run

# Overall function

final_output, ultra_rare, final_output_raw, ultra_rare_raw = load_variants(vcf_metrics_tsv, ultra_rare_variants_tsv, var_type)

final_output, ultra_rare, merged_output = filter_variants(final_output, ultra_rare, final_output_raw, ultra_rare_raw)
frac=min(1, prop_dn*final_output.shape[0]/ultra_rare.shape[0])
merged_output = pd.concat([final_output, ultra_rare.sample(frac=frac)])\
                .reset_index(drop=True)

if known_vars_exist:
    known_vars = pd.read_csv(known_vars_uri, sep='\t').set_index('VarKey', drop=False)
    merged_output['is_known'] = merged_output.VarKey.isin(known_vars.VarKey)

rep_reg_ur = pd.read_csv(ultra_rare_rep_regions, sep='\t', header=None)
rep_reg = pd.read_csv(rep_regions, sep='\t', header=None)

merged_output['repetitive_region'] = merged_output.VarKey.isin(pd.concat([rep_reg_ur, rep_reg])[3])

merged_output['multiallelic'] = (merged_output.DPC_sample!=merged_output.DP_sample)\
                |(merged_output.DPC_mother!=merged_output.DP_mother)\
                |(merged_output.DPC_father!=merged_output.DP_father)

# variant-level features
if len(variant_features) > 0:
    numeric = variant_features

    merged_output = merged_output[~merged_output[numeric].isna().any(axis=1)].reset_index(drop=True)

    X = merged_output[numeric].sample(frac=1)
    y = merged_output['label'].loc[X.index]

    results_variant_level, estimators_variant_level, importances_variant_level, oob_scores_variant_level = run_baggingPU_level_features(X, y, merged_output, numeric,
                                     n_estimators_rf, n_bags, metric, 'variant_level')
else:
    results_variant_level, importances_variant_level, oob_scores_variant_level = (pd.DataFrame() for _ in range(3))

# sample-level features
if len(sample_features) > 0:
    passes_variant_features = results_variant_level[(results_variant_level.pred_bag_optimized_variant_level==1)].VarKey
    merged_output = merged_output[merged_output.VarKey.isin(passes_variant_features)].reset_index(drop=True)

    numeric = sample_features

    merged_output = merged_output[~merged_output[numeric].isna().any(axis=1)].reset_index(drop=True)

    X = merged_output[numeric].sample(frac=1)
    y = merged_output['label'].loc[X.index]

    results_sample_level, estimators_sample_level, importances_sample_level, oob_scores_sample_level = run_baggingPU_level_features(X, y, merged_output, numeric,
                                     n_estimators_rf, n_bags, metric, 'sample_level')
else:
    results_sample_level, importances_sample_level, oob_scores_sample_level = (pd.DataFrame() for _ in range(3))

results_optimized = pd.concat([results_variant_level.set_index('VarKey', drop=False), 
           results_sample_level.set_index('VarKey', drop=False)[np.setdiff1d(results_sample_level.columns, results_variant_level.columns)]], axis=1)

feature_cols = [col for col in results_optimized.columns if 'feature' in col]
results_optimized['features'] = results_optimized[feature_cols].astype(str).agg(', '.join, axis=1)

importances_optimized = pd.concat([importances_variant_level, importances_sample_level], axis=1)

results_optimized.to_csv(f"{cohort_prefix}_{var_type}_{metric}_RF_results.tsv", sep='\t', index=False)
importances_optimized.to_csv(f"{cohort_prefix}_{var_type}_{metric}_RF_feature_importances.tsv", sep='\t', index=False)
# oob_scores_optimized.to_csv(f"{cohort_prefix}_{var_type}_{metric}_RF_oob_scores.tsv", sep='\t', index=False)
