# from wgs_bagging_pu_cross_validation.ipynb in GMKF-denovo-snvs-indels Analyses
import sys
import pandas as pd
import numpy as np
import sklearn.metrics
import sklearn.model_selection
import ast
import pickle 
import ast

from sklearn.ensemble import RandomForestClassifier
import xgboost as xgb
from baggingPU import BaggingClassifierPU

# https://heartbeat.comet.ml/positive-and-unlabelled-learning-recovering-labels-for-data-using-machine-learning-59c1def5452f
# https://roywrightme.wordpress.com/2017/11/16/positive-unlabeled-learning/

vcf_metrics_tsv = sys.argv[1]
ultra_rare_variants_tsv = sys.argv[2]
polyx_vcf = sys.argv[3]
cohort_prefix = sys.argv[4]
var_type = sys.argv[5]
numeric = sys.argv[6]
outlier_samples = sys.argv[7].split(',')
vqslod_cutoff = float(sys.argv[8])
model_type = sys.argv[9]  # string 'RF', 'XGB'
prop_dn = sys.argv[10]
base_params = pd.read_csv(sys.argv[11], header=None, index_col=0)[0].apply(ast.literal_eval).to_dict()  # dict
params = pd.read_csv(sys.argv[12], header=None, index_col=0)[0].apply(ast.literal_eval).to_dict()  # for cv
    
eval_metric = sklearn.metrics.f1_score

if model_type == 'RF':
    model_func = RandomForestClassifier    
if model_type == 'XGB':
    model_func = xgb.XGBClassifier

# Functions
def load_variants(vcf_metrics_tsv, ultra_rare_variants_tsv, polyx_vcf, var_type): 
    ultra_rare = pd.read_csv(ultra_rare_variants_tsv, sep='\t')
    ultra_rare = ultra_rare[ultra_rare.TYPE==var_type]
    # ultra_rare = ultra_rare[~ultra_rare.SAMPLE.isin(outlier_samples)]

    final_output = pd.read_csv(vcf_metrics_tsv, sep='\t')
    final_output = final_output[final_output.TYPE==var_type]
    # final_output = final_output[~final_output.SAMPLE.isin(outlier_samples)]

    polyx_df = pd.read_csv(polyx_vcf, sep='\t', comment='#', header=None)
    ultra_rare['POLYX'] = polyx_df[7].str.split('=').str[-1].astype(int)

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

    # TODO: remove when filter-rare-variants-hail is fixed?
    # ultra_rare = ultra_rare[ultra_rare.GQ_sample>=99]
    # ultra_rare = ultra_rare[(ultra_rare.GQ_mother>=30)&(ultra_rare.GQ_father>=30)]

    ultra_rare = ultra_rare[ultra_rare.POLYX<=10]

    try:
        ultra_rare = ultra_rare[~ultra_rare.VQSLOD.isna()]
        ultra_rare = ultra_rare[ultra_rare.VQSLOD>vqslod_cutoff]
    except:
        pass
    ultra_rare = ultra_rare[(~ultra_rare.VAF_sample.isna())&
                            (~ultra_rare.VAF_mother.isna())&
                            (~ultra_rare.VAF_father.isna())].reset_index(drop=True)

    ultra_rare_counts = pd.concat([ultra_rare_raw.groupby('SAMPLE').TYPE.value_counts(), 
               ultra_rare.groupby('SAMPLE').TYPE.value_counts()], axis=1)
    ultra_rare_counts.columns = ['before_filtering', 'after_filtering']
    ultra_rare_counts = pd.melt(ultra_rare_counts.reset_index(), id_vars=['SAMPLE', 'TYPE'], 
            value_vars=['before_filtering', 'after_filtering'],
            value_name='count', var_name='filter_status')

    ultra_rare['label'] = 1  # positive
    final_output['label'] = 0
    # overlapping variants set to unlabeled
    ultra_rare = ultra_rare[~ultra_rare.VarKey.isin(np.intersect1d(ultra_rare.VarKey, final_output.VarKey))]
  
    merged_output = pd.concat([ultra_rare, final_output]).reset_index(drop=True)
    merged_output['var_type'] = merged_output.label.map({1: 'ultra-rare', 0: 'unlabeled'})

    return final_output, ultra_rare, merged_output

def BaggingPU(X, y, kf, n_jobs=-1):
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
            n_estimators = 10, n_jobs = n_jobs, 
            max_samples = sum(y_train==0)  # all unlabeled
        )
        bc.fit(X_train, y_train)
        y_pred_bag[test_index] = bc.predict(X_test)
        output_bag[test_index] = bc.predict_proba(X_test)[:,1]
        classifiers.append(bc)
    return y_pred_bag, output_bag, classifiers

final_output, ultra_rare, final_output_raw, ultra_rare_raw = load_variants(vcf_metrics_tsv, ultra_rare_variants_tsv, polyx_vcf, var_type)

final_output, ultra_rare, merged_output = filter_variants(final_output, ultra_rare, final_output_raw, ultra_rare_raw)
merged_output = pd.concat([final_output, ultra_rare.sample(frac=min(1, prop_dn*final_output.shape[0]/ultra_rare.shape[0]))])\
                .reset_index(drop=True)

if numeric != 'false':
    numeric = numeric.split(',')
elif var_type == 'Indel':
    numeric = ['BaseQRankSum', 'MQ', 'MQRankSum', 'QD', 'GQ_parent'] 
elif var_type == 'SNV':
    numeric = ['BaseQRankSum', 'FS', 'MQ', 'MQRankSum', 'QD', 'SOR', 'GQ_parent']
numeric = np.intersect1d(numeric, np.intersect1d(ultra_rare.columns, final_output.columns)).tolist()

if model_type=='RF':  # can't handle NA features
    merged_output = merged_output[~merged_output[numeric].isna().any(axis=1)].reset_index(drop=True)

X = merged_output[numeric].sample(frac=1)
y = merged_output['label'].loc[X.index]
dtrain = xgb.DMatrix(data=X, label=y)

# Cross-validation
def my_score_func(predt: np.ndarray, dtrain: xgb.DMatrix):
    y = dtrain.get_label()
    y_pred = np.round(predt)
    FP = ((y==0) & (y_pred==1)).sum()
    TP = ((y==1) & (y_pred==1)).sum()
    FN = ((y==1) & (y_pred==0)).sum()
    TN = ((y==0) & (y_pred==0)).sum()
    FN_rate = FN / (FN + TN)
    FP_rate = FP / (FP + TP)
    if np.isnan(FP_rate):
        FP_rate = FN_rate
    return 'error', float(FN_rate / FP_rate)

def fit(alg, cv_folds=5, early_stopping_rounds=50):
    xgb_param = alg.get_xgb_params()
    xgb_cv = xgb.cv(xgb_param, dtrain, num_boost_round=alg.get_params()['n_estimators'], nfold=cv_folds,
            custom_metric=my_score_func, early_stopping_rounds=early_stopping_rounds, verbose_eval=False)
    return xgb_cv

def my_score_func_sk(y, y_pred):
    FP = ((y==0) & (y_pred==1)).sum()
    TP = ((y==1) & (y_pred==1)).sum()
    FN = ((y==1) & (y_pred==0)).sum()
    TN = ((y==0) & (y_pred==0)).sum()
    FN_rate = FN / (FN + TN)
    FP_rate = FP / (FP + TP)
    if np.isnan(FP_rate):
        FP_rate = FN_rate
    return float(FP_rate / FN_rate)

def cross_validation(params, recalc_n_estimators=3, early_stopping_rounds=50):
    # get n_estimators
    model = model_func(
        **base_params
    )
    
    new_params = base_params.copy()
    
    n_params = len(params)
    for i, (param, val) in enumerate(params.items()):
        param_msg = f"Optimizing parameter {i+1}/{n_params}, {param}..."
        print(param_msg)
        if (i%recalc_n_estimators==0):
            n_estimators = fit(model, early_stopping_rounds=early_stopping_rounds).shape[0]
            n_est_msg = f"Recalculated n_estimators to {n_estimators}"
            print(n_est_msg)
            new_params['n_estimators'] = n_estimators 
            
        model = model_func(**new_params)
        gsearch = sklearn.model_selection.GridSearchCV(estimator=model, param_grid={param: val},
                                scoring=sklearn.metrics.make_scorer(my_score_func_sk),
                                n_jobs=4, cv=5)
        gsearch.fit(X, y)
        print(gsearch.best_params_, gsearch.best_score_)
        new_params = new_params | gsearch.best_params_
        print('-'*len(param_msg))
        
    return model, new_params

if model_type=='RF':
    model, new_params = cross_validation(params, recalc_n_estimators=np.nan)
else:
    model, new_params = cross_validation(params)

finetuned_model = model_func(**new_params)
kf = sklearn.model_selection.KFold(n_splits=5)
y_pred_bag, output_bag, classifiers_bag = BaggingPU(X, y, kf, finetuned_model, n_jobs=-1)

print('---- {} ----'.format('Bagging PU'))
print(sklearn.metrics.confusion_matrix(y, y_pred_bag))

importances_bag = pd.DataFrame(np.array([[classifiers_bag[j].estimators_[i].feature_importances_ 
        for i in range(len(classifiers_bag[j].estimators_))] 
            for j in range(len(classifiers_bag))]).reshape((len(classifiers_bag)*len(classifiers_bag[0]),len(numeric))),
                columns=numeric)
    
results = pd.DataFrame({'label': y, 
                        'VarKey': merged_output.iloc[X.index].VarKey,
                        'predict_proba_bag': output_bag,
                        'pred_bag': y_pred_bag.astype(int)})

results.to_csv(f"{cohort_prefix}_baggingPU_{var_type}_{model_type}_CV_results.tsv", sep='\t', index=False)
importances_bag.to_csv(f"{cohort_prefix}_{var_type}_{model_type}_CV_feature_importances.tsv", sep='\t', index=False)
output = open(f"{cohort_prefix}_{var_type}_{model_type}_CV_model.pkl", 'wb')
pickle.dump(finetuned_model, output)
