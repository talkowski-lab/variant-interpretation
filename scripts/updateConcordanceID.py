#!/usr/bin/env python
#Author: Alba Sanchis Juan

"""
This script updates the ID of SVs that are concordant to the ID in the previous callset
"""

import argparse
import sys
import pysam

def main():
    """
    Parse arguments
    """
    parser = argparse.ArgumentParser(description='Parse arguments')
    parser.add_argument('--input', dest='vcfin', help='VCF input file')
    parser.add_argument('--output', dest='vcfout', help='VCF output file')
    args = parser.parse_args()

    if args.vcfin in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcfin)

    if args.vcfout in '- stdout'.split():
        out = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        out = pysam.VariantFile(args.vcfout, 'w', header=vcf.header)

    """
    Change ID column
    """
    for record in vcf:
        try:
            record.id = record.info["TRUTH_VID"]
            out.write(record)
        except:
            out.write(record)

if __name__ == '__main__':
    main()