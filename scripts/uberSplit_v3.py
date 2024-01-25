# Description: Split a large combined VCF into families as an input for trio de novo
# Usage: ulimit -n 1024; python uberSplit.py pedigree.txt VCF.txt.gz
# Author: Stephan Sanders
# Modified by Lindsay Liang 2020-4: Modify file handle part
# Modified by Shan Dong 2020-4: Add step to remove all homVar or "./." loci in output trio VCFs
# Modified by Shan Dong 2021-6-8: Recalcualte AC;AN;AF in output trio VCFs
# Modified by Michael Gilson 2021-7-26: Backport to Py2 to run in GCP pipeline. Can use 2to3 on this file to run in Python3. Will later use a py3 container in GCP. Fix str-int compare in phenotype parse.
# Modified by Shan Dong 2022-2: Dealing with pedigree with missing phenotype

#from pathlib import Path
#from rkstr8.clients.container import pipeline_task
#import csv
import gzip # open
import os
import sys # argv
#import subprocess
import math # ceil
from collections import defaultdict, namedtuple
# import re # 

class TrioStats:
    def __init__(self, hom_ref=0,missing=0,ok=0):
        self.hom_ref = hom_ref
        self.missing = missing
        self.ok = ok

pedFileName = sys.argv[1]
vcfFileName = sys.argv[2]
outPrefix = sys.argv[3]
statsFile = sys.argv[4]
BATCH_SIZE = float(sys.argv[5])

def ishomref(entry):
    gt = entry.split(':')[0]
    return gt == '0/0'
def ismissing(entry):
    gt = entry.split(':')[0]
    return gt == './.'

def get_ac_per_samp(entry):
    # entry looks like "0/1:21,19:40:99:431,0,510"
    gt = entry.split(":")[0]
    try:
        ac = int(gt.split("/")[0]) + int(gt.split("/")[1])
    except:
        ac = int(gt.split("|")[0]) + int(gt.split("|")[1])
    return ac

def update_ac(keyMetrics, chd_entry, pat_entry, mat_entry):
    fields = keyMetrics.split('\t')
    # based on VCF format INFO field is on 6th columm
    info_fields = fields[7].split(';')
    ac_new = get_ac_per_samp(chd_entry) + get_ac_per_samp(pat_entry) + get_ac_per_samp(mat_entry) 
    # AN will be 6 for trio, after removed  "./." genotype
    an_new = 6
    af_new = round(float(ac_new/an_new), 3)
    new_info = "AC=" + str(ac_new) + ";AF=" + str(af_new) + ";AN=" + str(an_new) + ";" + ";".join(info_fields[3:])
    return "\t".join(fields[:7]) + "\t" + new_info + "\t" + "\t".join(fields[8:]) # do not include newline

def write_trio_variant(fh, prefix, vcfCol, sampToCol, child, mother, father):
    # writes a VCF record for the trio's variant
    fh.write(
        '\t'.join((
            prefix, vcfCol[sampToCol[child]], vcfCol[sampToCol[father]], vcfCol[sampToCol[mother]]
        )) + '\n'
    )

def main():
    print('::RUNNING SS EXTRACT::')
    sys.stdout.flush()
    # Process the pedigree to get unique families
    sampToFam = {}
    famList = {}
    myTotalTrio = 0
    with open(pedFileName, 'rt') as inPed:
        for line in inPed:
            item = line.split('\t')
            famBasic = item[0] # the family id
            samp = item[1]
            father = item[2]
            mother = item[3]

            try:
                pheno = int(item[-1])
            except Exception as e:
                continue

            #
            # TODO: emit per-trio pedigrees here and can avoid the duplicate work I created.
            #

            # for each child
            if father != '0' and mother != '0':
                myTotalTrio += 1

                pheno = int(item[-1])
                if pheno == 2: # see https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format
                    #proband
                    fam = famBasic + '_trio_' + samp # fam is a uid for the trio (in a case where multiple trios are possible)
                elif pheno == 1:
                    #sibling
                    fam = famBasic + '_trio_' + samp
                elif pheno == 0:
                    # missing affected status
                    fam = famBasic + '_trio_' + samp
                else:
                    raise ValueError('Invalid phenotype {}'.format(pheno))
                famList[fam] = samp + '|' + father + '|' + mother

    with open(statsFile,'w') as stats_out:
        stats_out.write('## Splitter pedigree parsing\n')
        stats_out.write('#{}\n'.format('\t'.join(('fam_id','child','father','mother'))))
        for fam_id,trio_slug in famList.items():
            print(fam_id,trio_slug)
            stats_out.write(
                '{}\t{}\n'.format(
                    fam_id,'\t'.join(trio_slug.split('|'))
                )
            )

    safeFhCount = BATCH_SIZE # add trailing . to make this a float. Otherwise, because the division below for batchesNeeded 
    # is in an integer context, when myTotalTrio < safeFhCount, the result with be exactly 0, when we really want 1 batch.

    batchesNeeded =  int(math.ceil(myTotalTrio / safeFhCount)) # float is returned by math.ceil, invalid input to range
    #print('BN: ', str(batchesNeeded))

    #
    # Main. For each Batch..
    #
    # ensures output dir exists before attempting to write files within
    # if not os.path.exists(outPrefix):
    print(f"creating directory {outPrefix}")
    os.makedirs(outPrefix, exist_ok=True)   

    for batch in range(1, batchesNeeded + 1):
        print('Doing batch: '+str(batch)+' of '+str(batchesNeeded))
        # sys.stdout.flush()

        try:
            # let Python close the files for us, even in exceptional cases

            #
            # 1. Open a file for each trio
            #      

            famToFh = {}
            for famCount,fam in enumerate(famList):
                samp, father, mother = famList[fam].split('|')
                if famCount >= safeFhCount * (batch - 1) and famCount <= safeFhCount * batch:           
                    # make this dir/file entry
                    newFileName = '{}/{}.vcf'.format(outPrefix,fam)
                    famToFh[fam] = open(newFileName, 'w')

            with open(statsFile,'a') as stats_out:
                stats_out.write('## Splitter parsed pedigree to trio filenames(handles)\n')
                stats_out.write('#fam_id\tfile_name\n')
                for fam,fh in famToFh.items():
                    print(fam,fh.name)
                    stats_out.write('{}\t{}\n'.format(fam,fh.name))

            #
            # 2. Process the VCF to split into families
            #
            header = ''
            colOutList = []
            fhOutList = []

            filterStats = defaultdict(TrioStats)

            print("opening %s" % vcfFileName)

            with gzip.open(vcfFileName, 'rt') as in_vcf:
                for line in in_vcf:
                    if line.startswith('##'):
                        # VCF info, needs to go in each file
                        header = header + line
                    elif line.startswith('#CHROM'):    
                        # VCF column headings
                        vcfCol = line.rstrip().split('\t')
                        keyMetrics = '\t'.join(vcfCol[0:9])
                        header = header + keyMetrics

                        # map samples to columns number
                        colCount = 0
                        sampToCol = {}
                        for col in vcfCol:
                            sampToCol[col] = colCount
                            colCount += 1

                        # Print headers to each filehandle
                        for fam in famToFh:
                            fh = famToFh[fam]
                            child, father, mother = famList[fam].split('|')
                            write_trio_variant(fh, header, vcfCol, sampToCol, child, mother, father)

                    else:
                        # Rest of VCF, print out 
                        vcfCol = line.rstrip().split('\t')
                        keyMetrics = '\t'.join(vcfCol[0:9])

                        # For each trio (and its associated open filehandle for vcf data)
                        for fam in famToFh:
                            fh = famToFh[fam]
                            child, father, mother = famList[fam].split('|')
                            chd_entry = vcfCol[sampToCol[child]]
                            pat_entry = vcfCol[sampToCol[father]]
                            mat_entry = vcfCol[sampToCol[mother]]

                            if ishomref(chd_entry) and ishomref(pat_entry) and ishomref(mat_entry):
                                # update counter for stats
                                stats = filterStats[fam]
                                stats.hom_ref +=  1
                                continue
                            elif ismissing(chd_entry) or ismissing(pat_entry) or ismissing(mat_entry):
                                # update counter for stats
                                stats = filterStats[fam]
                                #print(stats.__dict__.keys())
                                stats.missing += 1
                                continue
                            else:
                                # update counter for stats
                                stats = filterStats[fam]
                                stats.ok += 1
                                new_keyMetrics = update_ac(keyMetrics,chd_entry,pat_entry,mat_entry)
                                write_trio_variant(fh, keyMetrics, vcfCol, sampToCol, child, mother, father)
        finally:
            # Clean up open file handles
            print('cleaning up')
            for fam in famToFh:
                famToFh[fam].close()
                
        with open("trio_vcf_files.txt", "w") as f:
            f.write("\n".join([fh.name for fh in famToFh.values()]))
        
        with open(statsFile,'a') as stats_out:
            stats_out.write('## Splitter trio filtering\n')
            stats_out.write('#{}\n'.format('\t'.join(('fam_id','hom_ref_filtered','missing_filtered','ok'))))
            for fam_id,stats in filterStats.items():
                stats_out.write(
                    '{}\t{}\n'.format(
                        fam_id,'\t'.join((str(stats.hom_ref),str(stats.missing),str(stats.ok)))
                    )
                )

        print(':::FINISHED FAM {}:::'.format(fam))
        sys.stdout.flush()

if __name__ == '__main__':
    main()


