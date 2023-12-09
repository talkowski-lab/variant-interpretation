#!/usr/bin/env python
#Author: Alba Sanchis Juan

"""
This script updates the ID of SVs that are concordant to the ID in the previous callset
"""

import argparse
import sys
import vcf

def main():
    """
    Parse arguments
    """
    parser = argparse.ArgumentParser(description='Parse arguments')
    parser.add_argument('--input', dest='vcfin', help='VCF input file')
    parser.add_argument('--output', dest='vcfout', help='VCF output file')
    args = parser.parse_args()

    vcfin = args.vcfin
    vcfout = args.vcfout

    if vcfin == '-':
        vcf_reader = vcf.Reader(sys.stdin)
    else:
        vcf_reader = vcf.Reader(filename=vcfin)

    if vcfout == None:
        vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
    else:
        vcf_writer = vcf.Writer(open(vcfout, 'w'), vcf_reader)

    """
    Change ID column
    """
    for record in vcf_reader:

        try:
            record.ID = record.INFO["TRUTH_VID"]
            vcf_writer.write_record(record)
        except:
            vcf_writer.write_record(record)

if __name__ == '__main__':
    main()