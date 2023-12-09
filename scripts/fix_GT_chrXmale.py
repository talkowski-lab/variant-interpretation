#!/usr/bin/env python

"""
This script converts 0/1 calls in the chrX non-PAR regions in males to 1/1
Author: Alba Sanchis Juan (asanchis-juan@mgh.harvard.edu)
"""

import argparse
import pysam
import sys
import pandas as pd
import numpy as np

def getMales(ped):
    #Get male ids
    male_ids = ped[ped["sex"] == 1]["subject_id"].to_numpy()
    return(male_ids)

def getIDs(ids):
    #Get SV ids to modify from file
    #Obtained running: bedtools intersect -a input.vcf.gz -b non_PAR.bed -f 0.5 | cut -f3 > ids.txt
    sv_ids = ids[0].to_numpy()
    return(sv_ids)

def fixChrX(vcf, out, males, ids):
    # Iterate over records
    for record in vcf:
        # If region of interest
        if record.id in ids:
            # Iterate over samples
            for male in males:
                if record.samples[male]['GT'] == (0,1) or record.samples[male]['GT'] == (1,0):
                    record.samples[male]['GT'] = (1,1)
        # Write record to file
        out.write(record)

def main():
    """
    Parse arguments
    """
    parser = argparse.ArgumentParser(description='Change het calls in chrX non-PAR regions in males to 1/1')
    parser.add_argument('--vcf', dest='vcf', help='VCF input file')
    parser.add_argument('--out', dest='out', help='Output file')
    parser.add_argument('--ped', dest='ped', help='Pedigree file')
    parser.add_argument('--ids', dest='ids', help='SV IDs in the chrX non-PAR region')

    args = parser.parse_args()

    """
    Read files
    """
    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    if args.out in '- stdout'.split():
        out = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        out = pysam.VariantFile(args.out, 'w', header=vcf.header)

    ped = pd.read_csv(args.ped, sep = '\t')
    ids = pd.read_csv(args.ids, sep='\t', header = None)

    males = getMales(ped)
    chrX_ids = getIDs(ids)

    fixChrX(vcf, out, males, chrX_ids)


if __name__ == '__main__':
    main()