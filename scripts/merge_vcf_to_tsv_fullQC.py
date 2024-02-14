#!/usr/bin/env python
"""
Script to merge each trio Triodenovo output VCF with the Triodenovo input VCF, then output the tsv format results
CHR, POS, REF, ALT and sample ID for Triodenovo output VCF are be used.
Other information (i.e. QUAL INFO, FORMAT) were extract from Triodenovo input VCF.
New columns (FAM, SAMPLE, TYPE (SNV/Indel), LEN, ID(KEY)) are generated
Usage:
python merge_vcf_to_tsv_fullQC.py -d triodenovo_out/ -i triodenovo_input/ -p pedigree.txt -o dnm.tsv

Shan Dong 2021-8-14 This script is for checking purpose
Shan Dong 2021-8-20 Make this script to output parsed sample level QC (from FORMAT) for both parents and child; 
Shan Dong 2021-8-23 write results into a single tsv; caculate GQ_mean; AB
Shan Dong 2021-9-2 deal with the situation when famID and sampleID contains "-"
Shan Dong 2021-9-9 fix issue with column oder and missing columns in INFO field - give a list of fixed columns (reason: KOR vcf contains missing INFO columns)
Shan Dong 2021-9-9 Add "case" and "control" status column (get the pedigree file)
Shan Dong 2022-1-24 dealing the situation with sampleID or familyID contains "." (complex naming system in MSSNG data)
"""
import gzip
import os, glob
import argparse
import numpy as np
import pandas as pd
from typing import DefaultDict


def main():
    # Parse arguments
    parser = create_arg_parser()
    args = parser.parse_args()
    
    # Collect trio list infor
    sample_info, triokey_to_samp = parse_pedigree(args.ped_path)
    
    # Collect triodenovo input keys (family ID and sample ID)
    res_vcf_paths = glob.glob(os.path.join(args.res_vcf_dir, '*.vcf.gz'))
    trio_variant_keys = DefaultDict(pd.DataFrame)  # key: familyID-sampleID; value: list of "chr_pos_ref_alt"
    for res_vcf_path in res_vcf_paths:
        trio_key = os.path.basename(res_vcf_path).replace('.denovos.vcf.gz', '').replace('_HP_VAF', '')
        with gzip.open(res_vcf_path, 'rt') as vcf_file:
            for line in vcf_file:
                if line.startswith('#CHROM'):
                    sample = triokey_to_samp[trio_key]
                    mother = sample_info[sample]['MotherID']
                    father = sample_info[sample]['FatherID']
                    sample_dict = {'DQ_sample': sample, 'DQ_father': father, 'DQ_mother': mother}
                    sample_order = {k: line.rstrip('\n').split('\t').index(s) for (k, s) in sample_dict.items()}
                    format_col = line.rstrip('\n').split('\t').index('FORMAT')
                if line.startswith('#'):
                    continue
                columns = line.rstrip('\n').split('\t')
                var_id = ':'.join([columns[0], columns[1], columns[3], columns[4]])
                try:
                    dq_idx = columns[format_col].split(':').index('DQ')
                except Exception as e:
                    print(trio_key)
                    print(columns[format_col])
                for k, idx in sample_order.items():
                    trio_variant_keys[trio_key].loc[var_id, k] = columns[idx].split(':')[dq_idx]
    
    # Make output 
    # abnd_columns = ['culprit', 'ExcessHet', 'InbreedingCoeff','MQ0']
    in_vcf_paths = glob.glob(os.path.join(args.in_vcf_dir, '*.vcf'))
    out_file = open(args.out_path, 'a')
    conter = 0
    for in_vcf_path in in_vcf_paths:
        trio_key = os.path.basename(in_vcf_path).replace('.vcf', '').replace('_HP_VAF', '')
        sampleID = triokey_to_samp[trio_key]
        print("Processing:" + sampleID)
        variant_df = parse_trio_vcf(in_vcf_path, sample_info, sampleID).set_index('ID', drop=False)
        result_df = variant_df.loc[trio_variant_keys[trio_key].index]
        result_df[trio_variant_keys[trio_key].columns] = trio_variant_keys[trio_key]
        if conter == 0:
            result_df.to_csv(out_file, sep = '\t', index = False, na_rep = '.')
        else:
            result_df.to_csv(out_file, sep = '\t', index = False, header = False, na_rep = '.')
        conter += 1
        print("Finished trio number:" + str(conter))
    out_file.close()
        
def create_arg_parser() -> argparse.ArgumentParser:
    """ Create an argument parser for this script and return it """
    # Create a top-level argument parser
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-d', '--res_dir', dest='res_vcf_dir', required=True, type=str,
                        help='Directory for Triodenovo ressult VCF files (per trio)')
    parser.add_argument('-i', '--in_dir', dest='in_vcf_dir', required=True, type=str,
                        help='Directory for Triodenovo input VCF files (per trio)')
    parser.add_argument('-p', '--ped_path', dest='ped_path', required=True, type=str,
                        help='path for standard pedigree file')
    parser.add_argument('-o', '--out_path', dest='out_path', required=False, type=str,
                        help='path for output table'
                        '(Default: merged_dnms_LF_triodenovo.tsv)', default='merged_dnms_LF_triodenovo.tsv')
    return parser

def parse_trio_vcf(vcf_path: str, sample_info: dict, sampleID: str, abnd_colnames: list = None) -> pd.DataFrame:
    """ Parse the VCF file and make a pandas.DataFrame object listing the annotated variants.
    :param vcf_path: The path of the VCF file
    :param abnd_colnames: The list of abandoned INFO names
                        (Warning: Unavailable column names will be ignored.)
    :return: The DataFrame object listing annotated variants
    """
    variant_df_rows = []
    variant_df_colnames = []
    id_list = []
    var_type = []
    var_len = []

    famID = sample_info[sampleID]['FamID']
    status = sample_info[sampleID]['Phenotype']
    fatherID = sample_info[sampleID]['FatherID']
    motherID = sample_info[sampleID]['MotherID']
    
    # Parse the VCF file
    with open(vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):  # The comments
                if line.startswith('#CHROM'):  # The header
                    variant_df_colnames = line[1:].rstrip('\n').split('\t')
            else:
                variant_df_row = line.rstrip('\n').split('\t')
                variant_df_rows.append(variant_df_row)
                id_list.append(":".join([variant_df_row[0],
                                        variant_df_row[1],
                                        variant_df_row[3],
                                        variant_df_row[4]]))
                # add variant type and length
                ref_len = len(variant_df_row[3])
                alt_len = len(variant_df_row[4])
                if ref_len == 1 and alt_len == 1:
                    var_type.append('SNV')
                    var_len.append(0)
                else:
                    var_type.append('Indel')
                    var_len.append(abs(ref_len - alt_len))
                
    vcf_df = pd.DataFrame(variant_df_rows, columns=variant_df_colnames)
    vcf_df['POS'] = pd.to_numeric(vcf_df['POS'])
    vcf_df['ID'] = np.array(id_list)
    vcf_df['TYPE'] = np.array(var_type)
    vcf_df['LEN'] = np.array(var_len)
    vcf_df['FAM'] = famID
    vcf_df['SAMPLE'] = sampleID
    vcf_df['STATUS'] = status
    vcf_df['VarKey'] = vcf_df['ID'] + ":" + vcf_df['SAMPLE']

    # Parse the INFO field
    fixed_info_columns = ['END','AC','AF','AN','BaseQRankSum','ClippingRankSum','DP','FS','MLEAC','MLEAF','MQ','MQRankSum','POLYX','QD','ReadPosRankSum','SOR','VQSLOD','cohort_AC','CSQ']
    info_strs = vcf_df['INFO'].to_numpy()
    info_dicts = list(map(parse_info_str, info_strs))
    info_df = pd.DataFrame(info_dicts)
    info_df = info_df.reindex(columns=fixed_info_columns)

    # Parse FORMAT field with trio info
    # child format columns
    format_return_columns = ['GT','AB','AD','DP','GQ','PGT','PID','PL','VAF']
    format_heads = vcf_df['FORMAT'].to_numpy()
    format_child_dict = list(map(parse_format_str, format_heads, vcf_df[sampleID].to_numpy()))
    format_child_df = pd.DataFrame(format_child_dict)
    format_child_df = format_child_df.reindex(columns=format_return_columns).add_suffix('_sample')
    # father format columns
    format_return_columns = ['GT','AB','AD','DP','GQ','VAF']
    format_heads = vcf_df['FORMAT'].to_numpy()
    format_father_dict = list(map(parse_format_str, format_heads, vcf_df[fatherID].to_numpy()))
    format_father_df = pd.DataFrame(format_father_dict)
    format_father_df = format_father_df.reindex(columns=format_return_columns).add_suffix('_father')
    # mother format columns
    format_heads = vcf_df['FORMAT'].to_numpy()
    format_mother_dict = list(map(parse_format_str, format_heads, vcf_df[motherID].to_numpy()))
    format_mother_df = pd.DataFrame(format_mother_dict)
    format_mother_df = format_mother_df.reindex(columns=format_return_columns).add_suffix('_mother')
    format_df = pd.concat([format_child_df, format_father_df, format_mother_df], axis='columns')
    format_df['GQ_mean'] = format_df[['GQ_father','GQ_mother','GQ_sample']].astype('int64').mean(axis='columns')

    # Concatenate those DataFrames
    variant_df = pd.concat([vcf_df.drop(columns=['INFO','FORMAT', sampleID, fatherID, motherID]), info_df, format_df], axis='columns')
    
    # Remove columns 
    if abnd_colnames is not None:
        variant_df.drop(columns=abnd_colnames, inplace=True, errors='ignore')

    #variant_df = variant_df.reindex(sorted(variant_df.columns), axis=1)
    return variant_df

def parse_info_str(info_str: str) -> dict:
    """ Parse the string in the INFO field of the VCF file and make a dictionary """
    
    info_dict = {}
    key_value_pairs = info_str.split(';')

    for key_value_pair in key_value_pairs:
        if '=' in key_value_pair:
            key, value = key_value_pair.split('=', 1)
            info_dict[key] = value
        else:
            # skip single tag in INFO - i.e. NEGATIVE_TRAIN_SITE
            continue

    return info_dict

def parse_format_str(format_head: str, format_str: str) -> dict:
    """ Parse the sample entries info using FORMAT string, make a dictionary (key contains sample role -- father, mother, child) 
        Example FORMAT:  GT:AB:AD:DP:FT:GQ:MIN_DP:MQ0:PGT:PID:PL:RGQ:SB
        Example FORMAT string: 0/1:0.520000:11,10:21:PASS:99:.:.:0|1:9837232_TTAA_T:381,0,516
        sample_role should be one of these: sample; father; mother 
    """
    format_dict = {}
    format_columns = format_head.split(":")
    format_content = format_str.split(":")
    for index, value in enumerate(format_content):
        if format_columns[index] == 'AB':
            continue
        if format_columns[index] == 'AD':
            try:
                format_dict['AB'] = float(value.split(',')[1]) / (float(value.split(',')[0]) + float(value.split(',')[1]))
            except:
                format_dict['AB'] = np.nan
        format_dict[format_columns[index]] = value
    
    return format_dict

def parse_pedigree(ped_path: str):
    sample_info = DefaultDict(dict) # key: sampleID; value: dict stored fam,father,mother,sex,status
    triokey_to_samp = DefaultDict(dict) # key: triokey; value:sample
    with open(ped_path) as ped_file:
        for line in ped_file:
            if not line.startswith('Fam'): # skip header line of pedigree file
                tabs = line.rstrip('\n').split('\t')
                # just extract full trio
                if tabs[2] != '0' and tabs[3] != '0': 
                    triokey_to_samp[tabs[0] + '_trio_' + tabs[1]] = tabs[1]
                    sample_info[tabs[1]]['FamID']= tabs[0] 
                    sample_info[tabs[1]]['FatherID']= tabs[2]
                    sample_info[tabs[1]]['MotherID']= tabs[3] 
                    if tabs[4] == '1':
                        sample_info[tabs[1]]['Sex'] = 'male'
                    elif tabs[4] == '2':
                        sample_info[tabs[1]]['Sex'] = 'female'
                    else:
                        sample_info[tabs[1]]['Sex'] = 'unknown'
                    if tabs[5] == '1':
                        sample_info[tabs[1]]['Phenotype'] = 'control'
                    elif tabs[5] == '2':
                        sample_info[tabs[1]]['Phenotype'] = 'case'
                    else:
                        sample_info[tabs[1]]['Phenotype'] = 'missing'
    return sample_info, triokey_to_samp

if __name__ == '__main__':
    main()