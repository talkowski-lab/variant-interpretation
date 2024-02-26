#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import sys
import os

final_output_uri = sys.argv[1]
final_output = pd.read_csv(final_output_uri, sep='\t')

final_output['Indel_type'] = final_output.apply(lambda x: 'Insertion' if (len(x.ALT) - len(x.REF)) > 0 else 'Deletion', axis=1)
final_output['START'] = final_output.apply(lambda x:  x.POS if x.Indel_type=='Insertion' else x.POS - x.LEN, axis=1)
final_output['END'] = final_output.apply(lambda x:  x.POS + x.LEN if x.Indel_type=='Insertion' else x.POS, axis=1)

if 'VarKey' not in final_output.columns:
    final_output['VarKey'] = final_output[['CHROM', 'POS', 'REF', 'ALT']].astype(str).agg(':'.join, axis=1)

new_filename = os.path.basename(final_output_uri).split('.tsv')[0] + '.bed'
final_output[['CHROM', 'START', 'END', 'VarKey', 'TYPE', 'SAMPLE']].to_csv(new_filename, sep='\t', index=False)

