# from wgs_positive_unlabeled_learning.ipynb in GMKF-denovo-snvs-indels Analyses
import sys
import pandas as pd
import numpy as np
import sklearn.metrics
import sklearn.model_selection
import ast

from sklearn.ensemble import RandomForestClassifier

from baggingPU import BaggingClassifierPU

# https://heartbeat.comet.ml/positive-and-unlabelled-learning-recovering-labels-for-data-using-machine-learning-59c1def5452f
# https://roywrightme.wordpress.com/2017/11/16/positive-unlabeled-learning/

vcf_metrics_tsv = sys.argv[1]
ultra_rare_variants_tsv = sys.argv[2]
cohort_prefix = sys.argv[3]
var_type = sys.argv[4]
numeric = sys.argv[5]

# Functions
def load_variants(vcf_metrics_tsv, ultra_rare_variants_tsv, var_type): 
    ultra_rare = pd.read_csv(ultra_rare_variants_tsv, sep='\t')
    ultra_rare = ultra_rare[ultra_rare.TYPE==var_type]

    final_output = pd.read_csv(vcf_metrics_tsv, sep='\t')
    final_output = final_output[final_output.TYPE==var_type]

    ultra_rare = ultra_rare[(~ultra_rare.AD_father.isna()) & (~ultra_rare.AD_mother.isna())].reset_index(drop=True)

    ultra_rare['VarKey'] = ultra_rare.ID + ':' + ultra_rare.SAMPLE

    
    ultra_rare['AB_mother'] = ultra_rare.AD_mother.apply(ast.literal_eval).apply(min) / ultra_rare.AD_mother.apply(ast.literal_eval).apply(sum)
    ultra_rare['AB_father'] = ultra_rare.AD_father.apply(ast.literal_eval).apply(min) / ultra_rare.AD_father.apply(ast.literal_eval).apply(sum)
    ultra_rare['AB_sample'] = ultra_rare.AD_sample.apply(ast.literal_eval).apply(min) / ultra_rare.AD_sample.apply(ast.literal_eval).apply(sum)

    ultra_rare = ultra_rare[~ultra_rare.PL_sample.isna()]
    ultra_rare[['PL_sample_0.0', 'PL_sample_0.1', 'PL_sample_1.1']] = ultra_rare.PL_sample.str.strip('][').str.split(', ', expand=True).astype(int)
    ultra_rare['GQ_mean'] = ultra_rare[['GQ_sample', 'GQ_mother', 'GQ_father']].mean(axis=1)

    # https://gatk.broadinstitute.org/hc/en-us/articles/360040507131-FilterVcf-Picard-
    # Allele balance is calculated for heterozygotes as the number of bases supporting the least-represented allele over the total number of base observations. Different heterozygote genotypes at the same locus are measured independently
    final_output['AB_mother'] = final_output.AD_mother.apply(ast.literal_eval).apply(min) / final_output.AD_mother.apply(ast.literal_eval).apply(sum)
    final_output['AB_father'] = final_output.AD_father.apply(ast.literal_eval).apply(min) / final_output.AD_father.apply(ast.literal_eval).apply(sum)
    final_output['AB_sample'] = final_output.AD_sample.apply(ast.literal_eval).apply(min) / final_output.AD_sample.apply(ast.literal_eval).apply(sum)

    final_output_raw = final_output.copy();
    ultra_rare_raw = ultra_rare.copy();
    return final_output, ultra_rare, final_output_raw, ultra_rare_raw

def filter_variants(final_output, ultra_rare, final_output_raw, ultra_rare_raw):
    final_output = final_output[final_output['VQSLOD']!='.'].reset_index(drop=True)

    ultra_rare = ultra_rare[(~ultra_rare.VQSLOD.isna())&
                            (~ultra_rare.VAF_sample.isna())&
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
        RandomForestClassifier(n_estimators = 1000, random_state=42), 
            n_estimators = 10, n_jobs = n_jobs, 
            max_samples = sum(y_train==0)  # all unlabeled
    #         max_samples = int(0.8*sum(y_train==0))  # 80% of unlabeled
        #     max_samples = min(ultra_rare_snv.shape[0], final_output_snv.shape[0])  
        )
        bc.fit(X_train, y_train)
        y_pred_bag[test_index] = bc.predict(X_test)
        output_bag[test_index] = bc.predict_proba(X_test)[:,1]
        classifiers.append(bc)
    return y_pred_bag, output_bag, classifiers

def run_PU_bagging(merged_output, numeric, n_splits=5):
    X = merged_output[merged_output.TYPE=='Indel'][numeric].sample(frac=1)
    y = merged_output[merged_output.TYPE=='Indel']['label'].loc[X.index]
    
    kf = sklearn.model_selection.KFold(n_splits=n_splits)
    y_pred_bag, output_bag, classifiers_bag = BaggingPU(X, y, kf, n_jobs=-1)

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
    return results, importances_bag

if numeric != '':
    numeric = numeric.split(',')
elif var_type == 'Indel':
    numeric = ['BaseQRankSum', 'MQ', 'MQRankSum', 'QD'] 
elif var_type == 'SNV':
    numeric = ['BaseQRankSum', 'FS', 'MQ', 'MQRankSum', 'QD', 'SOR']

final_output, ultra_rare, final_output_raw, ultra_rare_raw = load_variants(vcf_metrics_tsv, ultra_rare_variants_tsv, var_type)
final_output, ultra_rare, merged_output = filter_variants(final_output, ultra_rare, final_output_raw, ultra_rare_raw)
                                                         
results, importances_bag = run_PU_bagging(merged_output, numeric)

results.to_csv(f"{cohort_prefix}_baggingPU_{var_type}_results.tsv", sep='\t', index=False)
importances_bag.to_csv(f"{cohort_prefix}_{var_type}_feature_importances.tsv", sep='\t')