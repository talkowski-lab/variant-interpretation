# from wgs_positive_unlabeled_learning.ipynb in GMKF-denovo-snvs-indels Analyses
import sys
import pandas as pd
import numpy as np
import sklearn.metrics
import sklearn.model_selection

from sklearn.ensemble import RandomForestClassifier

from baggingPU import BaggingClassifierPU

# https://heartbeat.comet.ml/positive-and-unlabelled-learning-recovering-labels-for-data-using-machine-learning-59c1def5452f
# https://roywrightme.wordpress.com/2017/11/16/positive-unlabeled-learning/

vcf_metrics_tsv = sys.argv[1]
ultra_rare_variants_tsv = sys.argv[2]
cohort_prefix = sys.argv[3]
AC_threshold=int(sys.argv[4])
AF_threshold=float(sys.argv[5]) 
csq_AF_threshold=float(sys.argv[6])

# Functions
def load_variants(vcf_metrics_tsv, ultra_rare_variants_tsv): 
    ultra_rare = pd.read_csv(ultra_rare_variants_tsv, sep='\t')
    final_output = pd.read_csv(vcf_metrics_tsv, sep='\t')

    # ultra_rare['MLEAC'] = ultra_rare.MLEAC.str.strip('][').str.split(', ').str[0]
    # ultra_rare['MLEAF'] = ultra_rare.MLEAF.str.strip('][').str.split(', ').str[0]
    ultra_rare['CSQ'] = ultra_rare.CSQ.str.strip('][').str.split(', ').str[0]
    ultra_rare = ultra_rare[(~ultra_rare.AD_father.isna()) & (~ultra_rare.AD_mother.isna())].reset_index(drop=True)

    ultra_rare['VarKey'] = ultra_rare.ID + ':' + ultra_rare.SAMPLE

    ultra_rare['AB_mother'] = pd.Series([min(x.strip('][').split(', ')) for x in ultra_rare.AD_mother]).astype(int) / ultra_rare.DP_mother
    ultra_rare['AB_father'] = pd.Series([min(x.strip('][').split(', ')) for x in ultra_rare.AD_father]).astype(int) / ultra_rare.DP_father
    ultra_rare['AB_sample'] = pd.Series([min(x.strip('][').split(', ')) for x in ultra_rare.AD_sample]).astype(int) / ultra_rare.DP_sample

    ultra_rare = ultra_rare[~ultra_rare.PL_sample.isna()]
    ultra_rare[['PL_sample_0.0', 'PL_sample_0.1', 'PL_sample_1.1']] = ultra_rare.PL_sample.str.strip('][').str.split(', ', expand=True).astype(int)
    ultra_rare['GQ_mean'] = ultra_rare[['GQ_sample', 'GQ_mother', 'GQ_father']].mean(axis=1)
    final_output[['PL_sample_0.0', 'PL_sample_0.1', 'PL_sample_1.1']] = final_output.PL_sample.str.split(",", expand=True).astype(int)

    # np.setdiff1d(ultra_rare.columns, final_output.columns)
    # np.setdiff1d(final_output.columns, ultra_rare.columns)

    final_output = final_output[final_output.POLYX <= 10].reset_index(drop=True)

    # ultra_rare.TYPE.value_counts()

    # https://gatk.broadinstitute.org/hc/en-us/articles/360040507131-FilterVcf-Picard-
    # Allele balance is calculated for heterozygotes as the number of bases supporting the least-represented allele over the total number of base observations. Different heterozygote genotypes at the same locus are measured independently
    # I'm pretty sure merge_vcf_to_tsv_fullQC.py in step05 is calculating AB wrong (as VAF instead)...
    final_output['AB_mother'] = pd.Series([min(x.split(',')) for x in final_output.AD_mother]).astype(int) / final_output.DP_mother
    final_output['AB_father'] = pd.Series([min(x.split(',')) for x in final_output.AD_father]).astype(int) / final_output.DP_father
    final_output['AB_sample'] = pd.Series([min(x.split(',')) for x in final_output.AD_sample]).astype(int) / final_output.DP_sample

    # final_output.index = final_output.VarKey
    # ultra_rare.index = ultra_rare.VarKey

    final_output_raw = final_output.copy();
    ultra_rare_raw = ultra_rare.copy();
    return final_output, ultra_rare, final_output_raw, ultra_rare_raw

def filter_variants(final_output, ultra_rare, final_output_raw, ultra_rare_raw,
                    AC_threshold=3, AF_threshold=0.005, csq_AF_threshold=0.01):
    final_output = final_output[final_output['VQSLOD']!='.'].reset_index(drop=True)
    # ultra_rare = ultra_rare[(~ultra_rare.VQSLOD.isna())&
    #                         (~ultra_rare.VAF_sample.isna())].reset_index(drop=True)
    ultra_rare = ultra_rare[(~ultra_rare.VQSLOD.isna())&
                            (~ultra_rare.VAF_sample.isna())&
                            (~ultra_rare.VAF_mother.isna())&
                            (~ultra_rare.VAF_father.isna())].reset_index(drop=True)

    final_output = final_output[(final_output.cohort_AC<=AC_threshold) | (final_output.AF<=AF_threshold)].reset_index(drop=True)
    final_output = final_output[final_output['PL_sample_0.1']==0].reset_index(drop=True)

    csq_columns = ["Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature",
                   "BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position",
                   "Protein_position","Amino_acids","Codons","Existing_variation","ALLELE_NUM",
                   "DISTANCE","STRAND","FLAGS","VARIANT_CLASS","MINIMISED","SYMBOL_SOURCE",
                   "HGNC_ID","CANONICAL","MANE_SELECT","MANE_PLUS_CLINICAL","TSL","APPRIS",
                   "CCDS","ENSP","SWISSPROT","TREMBL","UNIPARC","UNIPROT_ISOFORM","GENE_PHENO",
                   "SIFT","PolyPhen","DOMAINS","miRNA","HGVS_OFFSET","AF","AFR_AF","AMR_AF",
                   "EAS_AF","EUR_AF","SAS_AF","gnomADe_AF","gnomADe_AFR_AF","gnomADe_AMR_AF",
                   "gnomADe_ASJ_AF","gnomADe_EAS_AF","gnomADe_FIN_AF","gnomADe_NFE_AF","gnomADe_OTH_AF",
                   "gnomADe_SAS_AF","gnomADg_AF","gnomADg_AFR_AF","gnomADg_AMI_AF","gnomADg_AMR_AF",
                   "gnomADg_ASJ_AF","gnomADg_EAS_AF","gnomADg_FIN_AF","gnomADg_MID_AF","gnomADg_NFE_AF",
                   "gnomADg_OTH_AF","gnomADg_SAS_AF","MAX_AF","MAX_AF_POPS","CLIN_SIG","SOMATIC",
                   "PHENO","PUBMED","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE",
                   "TRANSCRIPTION_FACTORS","LoF","LoF_filter","LoF_flags","LoF_info"]

    csq = final_output.CSQ.str.split(',', expand=True)[0].str.split('|', expand=True)
    csq.index = final_output.index
    csq.columns = csq_columns

    final_output[['gnomADe_AF', 'MAX_AF', '1KG_AF']] = csq[['gnomADe_AF', 'MAX_AF', 'AF']]

    csq_ultra_rare = ultra_rare.CSQ.str.split(',', expand=True)[0].str.split('|', expand=True)
    csq_ultra_rare.index = ultra_rare.index
    csq_ultra_rare.columns = csq_columns
    ultra_rare[['gnomADe_AF', 'MAX_AF', '1KG_AF']] = csq_ultra_rare[['gnomADe_AF', 'MAX_AF', 'AF']]

    af_csq = csq[csq.MAX_AF!='']
    final_output = final_output.drop(af_csq[(af_csq.MAX_AF.astype(float)>=csq_AF_threshold)].index)

    final_output_counts = pd.concat([final_output_raw.groupby('SAMPLE').TYPE.value_counts(), 
               final_output.groupby('SAMPLE').TYPE.value_counts()], axis=1)
    final_output_counts.columns = ['before_filtering', 'after_filtering']
    final_output_counts = pd.melt(final_output_counts.reset_index(), id_vars=['SAMPLE', 'TYPE'], 
            value_vars=['before_filtering', 'after_filtering'],
            value_name='count', var_name='filter_status')

    # final_output_counts.groupby(['TYPE', 'filter_status'])['count'].median()

    ultra_rare_counts = pd.concat([ultra_rare_raw.groupby('SAMPLE').TYPE.value_counts(), 
               ultra_rare.groupby('SAMPLE').TYPE.value_counts()], axis=1)
    ultra_rare_counts.columns = ['before_filtering', 'after_filtering']
    ultra_rare_counts = pd.melt(ultra_rare_counts.reset_index(), id_vars=['SAMPLE', 'TYPE'], 
            value_vars=['before_filtering', 'after_filtering'],
            value_name='count', var_name='filter_status')

    # ultra_rare.TYPE.value_counts()

    ultra_rare['label'] = 1  # positive
    final_output['label'] = 0

    merged_output = pd.concat([ultra_rare, final_output]).reset_index(drop=True)
    merged_output['var_type'] = merged_output.label.map({1: 'ultra-rare', 0: 'unlabeled'})
    # merged_output.groupby('var_type').TYPE.value_counts().unstack()
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

def run_PU_bagging_indels(merged_output, numeric_indel, n_splits=5):
    X_indel = merged_output[merged_output.TYPE=='Indel'][numeric_indel].sample(frac=1)
    y_indel = merged_output[merged_output.TYPE=='Indel']['label'].loc[X_indel.index]
    
    kf = sklearn.model_selection.KFold(n_splits=n_splits)
    y_pred_bag_indel, output_bag_indel, classifiers_bag_indel = BaggingPU(X_indel, y_indel, kf, n_jobs=-1)

    print('---- {} ----'.format('Bagging PU'))
    print(sklearn.metrics.confusion_matrix(y_indel, y_pred_bag_indel))

    importances_bag_indel = pd.DataFrame(np.array([[classifiers_bag_indel[j].estimators_[i].feature_importances_ 
        for i in range(len(classifiers_bag_indel[j].estimators_))] 
            for j in range(len(classifiers_bag_indel))]).reshape((len(classifiers_bag_indel)*len(classifiers_bag_indel[0]),len(numeric_indel))),
                columns=numeric_indel)
    
    results_indel = pd.DataFrame({'label': y_indel, 
                            'VarKey': merged_output.iloc[X_indel.index].VarKey,
                            'pred_bag': y_pred_bag_indel.astype(int)})
    return results_indel, importances_bag_indel

numeric_indel = ['BaseQRankSum', 'MQ', 'MQRankSum', 'QD'] 
#                  'ReadPosRankSum', 'DP', 'SOR']
 
final_output, ultra_rare, final_output_raw, ultra_rare_raw = load_variants(vcf_metrics_tsv, ultra_rare_variants_tsv)
final_output, ultra_rare, merged_output = filter_variants(final_output, ultra_rare, final_output_raw, ultra_rare_raw,
                                                          AC_threshold=AC_threshold, AF_threshold=AF_threshold, csq_AF_threshold=csq_AF_threshold)
results_indel, importances_bag_indel = run_PU_bagging_indels(merged_output, numeric_indel)

results_indel.to_csv(f"{cohort_prefix}_baggingPU_results.tsv", sep='\t', index=False)
importances_bag_indel.to_csv(f"{cohort_prefix}_feature_importances.tsv", sep='\t')