#!/usr/bin/env python

"""
Author: Alba Sanchis-Juan
Description: SV de novo filtering script
"""

# Import libraries
import argparse
import yaml
import math
import numpy as np
import pandas as pd
import pybedtools
pd.options.mode.chained_assignment = None  # default='warn'

# Define functions
def verbosePrint(msg, verbose):
    if verbose == "True":
        print(msg)

def getCount(row, list):
    return (
        sum(samp in list for samp in row['samples'].split(','))
    )

def getParentsFrequency(row, list):
    return (
        sum(samp in list for samp in row['samples'].split(','))/len(list)
    )

def getFamilyCount(row, ped):
    sample = row['sample']
    father = ped[(ped['subject_id'] == sample)]['paternal_id'].values[0]
    mother = ped[(ped['subject_id'] == sample)]['maternal_id'].values[0]
    parents = [father, mother]
    return (
        sum(samp in parents for samp in row['samples'].split(','))
    )

def variantInfo(row, field, vcf):
    row_name = row['name']
    samp = row['sample']

    filt_pos = vcf[(vcf['ID'] == row_name)]['FORMAT'].str.split(':').tolist()

    if field in filt_pos[0]:
        idx = filt_pos[0].index(field)
        if samp in vcf.columns:
            samp_info = vcf[(vcf['ID'] == row_name)][samp].str.split(':').tolist()
            samp_info_field = samp_info[0][idx]
        else:
            samp_info_field = 'NA'
    else:
        samp_info_field = 'NA'

    return(samp_info_field)

def addFamily(row, ped):
    sample = row['sample']
    fam = ped[(ped['subject_id'] == sample)]['family_id'].values
    chr = row['chrom']
    if fam.size != 0:
        return ('_'.join([str(fam[0]), chr]))


def getInfoParent(row, ped, vcf, parent, field):    
    row_name = row['name']
    sample = row['sample']
    parent_id = ped[(ped['subject_id'] == sample)][parent].values

    filt_pos = vcf[(vcf['ID'] == row_name)]['FORMAT'].str.split(':').tolist()

    if (field in filt_pos[0]):
        idx = filt_pos[0].index(field)
        parent_info = vcf[(vcf['ID'] == row_name)][parent_id[0]].str.split(':').tolist()
        parent_info_field = parent_info[0][idx]
    else:
        return('NA')

    return(parent_info_field)

def getParents(proband,ped):
    mother = ped[(ped['subject_id'] == proband)]['maternal_id'].values
    father = ped[(ped['subject_id'] == proband)]['paternal_id'].values
    parent_list = [mother,father]
    return(parent_list)

def getFamilyID(row,ped):
    samp = row['sample']
    family_id = ped[(ped['subject_id'] == samp)]['family_id'].values
    return(family_id)

    

def main():
    """
    Parse arguments
    """
    parser = argparse.ArgumentParser(description='Parse arguments')
    parser.add_argument('--bed', dest='bed', help='Input BED file')
    parser.add_argument('--ped', dest='ped', help='Ped file')
    parser.add_argument('--vcf', dest='vcf', help='VCF file')
    parser.add_argument('--disorder', dest='disorder', help='Genomic disorder regions')
    parser.add_argument('--out', dest='out', help='Output file')
    parser.add_argument('--raw_proband', dest='raw_proband', help='Directory with raw SV calls - output from m01')
    parser.add_argument('--raw_parents', dest='raw_parents', help='Directory with raw SV calls - output from m01')
    parser.add_argument('--outliers', dest='outliers', help='Output file with SV calls in outlier samples')
    parser.add_argument('--config', dest='config', help='Config file')
    parser.add_argument('--verbose', dest='verbose', help='Verbosity')
    args = parser.parse_args()

    bed_file = args.bed
    ped_file = args.ped
    vcf_file = args.vcf
    disorder_file = args.disorder
    out_file = args.out
    raw_file_proband = args.raw_proband
    raw_file_parent = args.raw_parents
    outlier_samples = args.outliers
    verbose = args.verbose
    config_file = args.config



    with open(config_file, "r") as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    large_cnv_size = int(config['large_cnv_size'])
    gnomad_col = config['gnomad_col']
    gnomad_AF = float(config['gnomad_AF'])
    parents_AF = float(config['parents_AF'])
    parents_SC = int(config['parents_SC'])
    large_raw_overlap = float(config['large_raw_overlap'])
    small_raw_overlap = float(config['small_raw_overlap'])
    cohort_AF = float(config['cohort_AF'])

    # Read files
    verbosePrint('Reading Input Files', verbose)
    bed = pd.read_csv(bed_file, sep='\t').replace(np.nan, '', regex=True)
    bed = bed[(bed['samples'] != "")]
    bed.rename(columns={'#chrom': 'chrom'}, inplace=True)
    #vcf = pd.read_csv(vcf_file, sep='\t', comment='#', names=('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'))
    vcf = pd.read_csv(vcf_file, sep='\t')
    ped = pd.read_csv(ped_file, sep='\t')
    disorder = pd.read_csv(disorder_file, sep='\t', header=None)
    raw_bed_colnames = colnames=['ID', 'start', 'end', 'svtype', 'sample'] 
    raw_bed_child = pd.read_csv(raw_file_proband, sep='\t', names= raw_bed_colnames, header=None).replace(np.nan, '', regex=True)
    raw_bed_parent = pd.read_csv(raw_file_parent, sep='\t', names= raw_bed_colnames, header=None).replace(np.nan, '', regex=True)



    # Get parents and children ids
    verbosePrint('Getting parents/children/affected/unaffected IDs', verbose)
    # families = ped[(ped['family_size'] >= 3) &
    families = ped[(ped['family_size'] == 3) &
                   (ped['paternal_id'] != "0") &
                   (ped['maternal_id'] != "0")]['family_id'].values
    trios = ped[(ped['family_id'].isin(families))]
    parents = trios[(trios['paternal_id'] == "0") & (trios['maternal_id'] == "0")]['subject_id'].values
    children = trios[(trios['paternal_id'] != "0") & (trios['maternal_id'] != "0")]['subject_id'].values

    # Get counts at cohort level and number of parents/children
    verbosePrint('Getting counts', verbose)
    bed['num_children'] = bed.apply(lambda r: getCount(r, children), axis=1)
    bed['num_parents'] = bed.apply(lambda r: getCount(r, parents), axis=1)
    bed['AF_parents'] = bed.apply(lambda r: getParentsFrequency(r, parents), axis=1)

    # Remove mCNVs, BNDs and SVs in sex chromosomes
    verbosePrint('Remove BND and mCNV', verbose)
    bed = bed[(~bed['svtype'].isin(['BND', 'CNV']))]
    # bed = bed[(~bed['svtype'].isin(['BND', 'CNV']) & (~bed['chrom'].isin(["chrY", "chrX"])))]

    
    # Flag if small or large CNV based on 5Kb cutoff
    verbosePrint('Flagging calls depending on size', verbose)
    bed['is_large_cnv'] = (bed['SVLEN'] >= large_cnv_size) & ((bed['svtype'] == 'DEL') | (bed['svtype'] == 'DUP'))
    bed['is_small_cnv'] = (bed['SVLEN'] < large_cnv_size) & ((bed['svtype'] == 'DEL') | (bed['svtype'] == 'DUP'))
    bed['is_depth_only'] = (bed['EVIDENCE'] == "RD")

    # Split into one row per sample
    verbosePrint('Split into one row per sample', verbose)
    bed_split = bed.assign(sample=bed.samples.str.split(",")).explode('sample')
   


    # Sepparate variants in children and parents
    verbosePrint('Sepparate variants in children and parents', verbose)
    bed_child = bed_split[bed_split['sample'].str.contains("|".join(children))]
    bed_parents = bed_split[bed_split['sample'].str.contains("|".join(parents))]


    print("Size of bed_child:",str(len(bed_child)))


    # Flag variants in genomic disorder (GD) regions
    verbosePrint('Flagging variants in GD', verbose)
    bed_child['in_gd'] = bed_child['name'].isin(disorder[0])

    # Filter out by frequency - AF gnomad < 0.01 OR inGD
    verbosePrint('Filtering by frequency', verbose)
    bed_child[gnomad_col] = pd.to_numeric(bed_child[gnomad_col])
    bed_child["AF"] = pd.to_numeric(bed_child["AF"])

    bed_child = bed_child[( ((bed_child[gnomad_col] <= gnomad_AF) & (bed_child['AF'] <= cohort_AF)) |
                           ((bed_child[gnomad_col].isnull() )  & (bed_child['AF'] <= cohort_AF)) |
                          (bed_child['in_gd'] == True) )]  


    

    print("Size of bed_child after filtering by frequency:",str(len(bed_child)))

    # Get counts within family and remove if SV in parents
    verbosePrint('Keep variants in children only', verbose)
    bed_child['num_parents_family'] = bed_child.apply(lambda r: getFamilyCount(r, ped), axis=1)
    bed_child = bed_child[ (bed_child['num_parents_family'] == 0) &
                           (bed_child['num_children'] >= 1) &
                           (bed_child['AF_parents'] <= parents_AF) &
                           (bed_child['num_parents'] <= parents_SC) ]

    print("Size of bed_child after keeping variants in children only:",str(len(bed_child)))



    #print(bed_child['svtype'].value_counts()['INS']) #JUMPED DOWN TO 14

    # Extract info from the VCF file - no PE_GT, PE_GQ, SR_GT, SR_GQ
    verbosePrint('Appending FILTER information', verbose)
    bed_child['GT']    = bed_child.apply(lambda r: variantInfo(r, 'GT', vcf), axis=1)
    bed_child['EV']    = bed_child.apply(lambda r: variantInfo(r, 'EV', vcf), axis=1)
    bed_child['GQ']    = bed_child.apply(lambda r: variantInfo(r, 'GQ', vcf), axis=1)
    bed_child['RD_CN'] = bed_child.apply(lambda r: variantInfo(r, 'RD_CN', vcf), axis=1)
    bed_child['RD_GQ'] = bed_child.apply(lambda r: variantInfo(r, 'RD_GQ', vcf), axis=1)
    bed_child['PE_GQ'] = bed_child.apply(lambda r: variantInfo(r, 'PE_GQ', vcf), axis=1)
    bed_child['PE_GT'] = bed_child.apply(lambda r: variantInfo(r, 'PE_GT', vcf), axis=1)
    bed_child['SR_GQ'] = bed_child.apply(lambda r: variantInfo(r, 'SR_GQ', vcf), axis=1)
    bed_child['SR_GT'] = bed_child.apply(lambda r: variantInfo(r, 'SR_GT', vcf), axis=1)


    # Remove WHAM only and GT = 1
    verbosePrint('Remove wham only and GT=1 calls', verbose)
    bed_child = bed_child[~((bed_child['ALGORITHMS'] == "wham") & (bed_child['GQ'] == '1'))]

    print("Size of bed_child after removing wham only and GT=1 calls:",str(len(bed_child)))


    # LARGE CNV: Check for false negative in parents: check depth in parents, independently of the calls
    verbosePrint('Large CNVs check', verbose)
    # 1. Add parental RD_CN field and strip out if same as in proband
    bed_child['paternal_rdcn'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'paternal_id', 'RD_CN'), axis=1)
    bed_child['maternal_rdcn'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'maternal_id', 'RD_CN'), axis=1)
    bed_child['paternal_srgq'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'paternal_id', 'SR_GQ'), axis=1)
    bed_child['maternal_srgq'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'maternal_id', 'SR_GQ'), axis=1)

    bed_child = bed_child.loc[(bed_child['is_large_cnv'] == False) | ((bed_child['RD_CN'] != bed_child['maternal_rdcn']) & (bed_child['RD_CN'] != bed_child['paternal_rdcn']))]
    
    print("Size of bed_child after removing if RD_CN field is same as parent:",str(len(bed_child)))

    # 2. Check if call in parents with bedtools coverage (332)
    verbosePrint('CNV present in parents check', verbose)
    cols_keep = ['family_chrom', 'start', 'end', 'name', 'svtype', 'sample']
    bed_child_large = bed_child[(bed_child['is_large_cnv'] == True)]
    bed_child_large['family_chrom'] = bed_child_large.apply(lambda r: addFamily(r, ped), axis=1)
    bed_child_large = bed_child_large[cols_keep].to_string(header=False, index=False)
    bed_child_large = pybedtools.BedTool(bed_child_large, from_string=True).sort()

    bed_parents_large = bed_parents[(bed_parents['is_large_cnv'] == True)]
    bed_parents_large['family_chrom'] = bed_parents_large.apply(lambda r: addFamily(r, ped), axis=1)
    bed_parents_large = bed_parents_large[cols_keep].to_string(header=False, index=False)
    bed_parents_large = pybedtools.BedTool(bed_parents_large, from_string=True).sort()

    bed_overlap = bed_child_large.coverage(bed_parents_large).to_dataframe(disable_auto_names=True, header=None)
    names_overlap = bed_overlap[(bed_overlap[9] >= 0.5)][3].to_list()

    ##Add if overlap to bed child table
    bed_child['overlap_parent'] = (bed_child['name'].isin(names_overlap))



    # Small calls:
    # If RD,SR and <1Kb, treat RD,SR as SR
    verbosePrint('Small CNVs check', verbose)
    bed_child['EVIDENCE_FIX'] = bed_child['EVIDENCE']
    bed_child[(bed_child['SVLEN'] <= large_cnv_size) & \
                    (bed_child['EVIDENCE'] == "RD,SR") & \
                    ((bed_child['svtype'] == 'DEL') | (bed_child['svtype'] == 'DUP'))]['EVIDENCE_FIX'] = "SR"

    ##Add if variant contains RD evidence
    bed_child['contains_RD'] = bed_child.EVIDENCE_FIX.str.contains('RD')


    ###############
    ## FILTERING ##
    ###############
    verbosePrint('Filtering out calls', verbose)
    # Keep if in GD region
    keep_gd = bed_child[(bed_child['in_gd'] == True)]['name'].to_list()
    print(len(keep_gd))
    # Filter out if large CNVs don't have RD support and parents overlap
    keep_large = bed_child[ (bed_child['is_large_cnv'] == True) &
                            (bed_child['contains_RD'] == True) &
                            (bed_child['overlap_parent'] == False) ]['name'].to_list()
    print("large_cnvs")
    print(len(keep_large))
    # Filter out if small cnvs that are SR-only don't have BOTHSIDES_SUPPORT
    keep_small = bed_child[(bed_child['is_small_cnv'] == True) &
                           ( (bed_child['EVIDENCE_FIX'] != 'SR') |
                           ( (bed_child['EVIDENCE_FIX'] == 'SR') &
                             (bed_child.FILTER.str.contains('BOTHSIDES_SUPPORT'))))
                          ]['name'].to_list()
    print("small_cnvs before filtering:")
    print(bed_child[(bed_child['is_small_cnv'] == True)]["FILTER"])
    print(bed_child[(bed_child['is_small_cnv'] == True)])
    print("small cnvs after filtering")
    print(keep_small)

    # Keep any other SV type
    keep_other_sv = bed_child[~bed_child['SVTYPE'].isin(['DEL', 'DUP'])]['name'].to_list()

    # Join lists and keep unique values
    keep_names = list(set(keep_gd + keep_large + keep_small + keep_other_sv))

    #Subset table
    bed_filt = bed_child[(bed_child['name'].isin(keep_names))]

    print("Size of bed_filt after removing large CNVs that dont have RD support and parents overlap and small cnvs that are SR only and dont have support on both sides:",str(len(bed_filt)))
    

    bed_filt.to_csv(path_or_buf="bed_filt.bed", mode='a', index=False, sep='\t', header=True)
    

    

    #####################
    ## REMOVE OUTLIERS ##
    #####################
    verbosePrint('Removing outliers', verbose)
    #1.5*IQR + 3Q
    sample_counts = bed_filt['sample'].value_counts()
    #print(sample_counts)
    q3, q1 = np.percentile(sample_counts.tolist(), [75, 25])
    iqr = q3 - q1
    count_threshold = int(math.ceil((1.5*iqr + q3)*3))
    #print(count_threshold) #multiplying by 3 since threshold seemed too strict
    #exit()
    sample_counts2 = sample_counts.rename_axis('sample').reset_index(name='count')
    samples_keep = sample_counts2[(sample_counts2['count'] <= count_threshold)]['sample'].tolist()

    #print(samples_keep)
    ###################
    ##RAW FILES CHECK##
    ###################
    #Keep all CPT and CTX - they are rare and the type is different than in the raw files
    #Check only DEL,DUP and INS for raw evidence
    verbosePrint('Checking raw files', verbose)

    #Reformat raw files
    #print(raw_bed)
    raw_bed_ref_child = raw_bed_child.to_string(header=False, index=False)
    raw_bed_ref_child = pybedtools.BedTool(raw_bed_ref_child, from_string=True)


    raw_bed_ref_parent = raw_bed_parent.to_string(header=False, index=False)
    raw_bed_ref_parent = pybedtools.BedTool(raw_bed_ref_parent, from_string=True)



    #raw_bed_ref = pybedtools.BedTool.from_dataframe(raw_bed, disable_auto_names=True)
    #print(raw_bed_ref)
    

    #raw_bed_ref_child = pybedtools.BedTool(raw_file_proband)
    #raw_bed_ref_parent = pybedtools.BedTool(raw_file_parent)
    
    #name = svtype
    #score = sample
    #print(raw_bed_ref[1].score)
    #print(raw_bed_ref.filter(lambda x: x.score == '__3176003__2f003a'))
    #it has the mom and child in raw_file
    #exit()


    #Reformat de novo filt calls

    bed_filt['family_id'] = bed_filt.apply(lambda r: getFamilyID(r,ped), axis=1)
    bed_filt['chrom_type_sample'] = bed_filt['chrom'] + "_" + bed_filt['SVTYPE'] + "_" + bed_filt['sample']
    bed_filt['chrom_type_family'] = bed_filt['chrom'] + "_" + bed_filt['SVTYPE'] + "_" + bed_filt['family_id'].astype(str).str.strip("[]")



    cols_keep2 = ['chrom_type_sample', 'start', 'end', 'name', 'svtype', 'sample']
    cols_keep3 = ['chrom_type_family', 'start', 'end', 'name', 'svtype', 'sample']


    ##INSERTIONS: dist < 100bp (& INS ratio > 0.1 & INS ratio < 10 ?) (dis is START-start and INS ratio is SVLEN/svlen)
    verbosePrint('Checking insertions in raw files', verbose)

    verbosePrint('Checking if insertions in proband is in raw files', verbose)

    bed_filt_ins = bed_filt[bed_filt['SVTYPE'] == 'INS']

    print("Number of insertions found:",str(len(bed_filt_ins)))

    if (len(bed_filt_ins.index) > 0):
        bed_filt_ins_proband = bed_filt_ins[cols_keep2].to_string(header=False, index=False)
        bed_filt_ins_proband = pybedtools.BedTool(bed_filt_ins_proband, from_string=True).sort()
        # bed_filt_ins_ref = pybedtools.BedTool(bed_filt_ins_ref, from_string=True)
        bed_filt_ins_overlap_proband = bed_filt_ins_proband.closest(raw_bed_ref_child).to_dataframe(disable_auto_names=True, header=None)
        bed_filt_ins_overlap_proband['is_close'] = abs(bed_filt_ins_overlap_proband[7] - bed_filt_ins_overlap_proband[1]) < 100

        ins_names_overlap_proband = bed_filt_ins_overlap_proband[(bed_filt_ins_overlap_proband['is_close'] == True)][3].to_list()
    else:
        ins_names_overlap_proband = ['']
    print("Number of insertions supported by raw evidence:",str(len(ins_names_overlap_proband)))
    print(ins_names_overlap_proband)

    #subset_parents = a.filter(lambda x: x.name == )
    #We will want to compare original raw evidence with subsetted raw evidence and check if parent has variant
    #print(bed_filt_ins_overlap_proband)
    #print(ins_names_overlap_proband)




    verbosePrint('Checking if insertions in proband are also in raw files for the parents', verbose)

    if (len(bed_filt_ins.index) > 0):
        bed_filt_ins_fam = bed_filt_ins[ cols_keep3 ].to_string(header=False, index=False)
        bed_filt_ins_fam = pybedtools.BedTool(bed_filt_ins_fam, from_string=True).sort()
        # bed_filt_ins_ref = pybedtools.BedTool(bed_filt_ins_ref, from_string=True)
        bed_filt_ins_overlap_parents = bed_filt_ins_fam.closest(raw_bed_ref_parent).to_dataframe(disable_auto_names=True, header=None)
        bed_filt_ins_overlap_parents[ 'is_close' ] = abs(bed_filt_ins_overlap_parents[ 7 ] - bed_filt_ins_overlap_parents[ 1 ]) < 100
        

        ins_names_overlap_parents = bed_filt_ins_overlap_parents[ (bed_filt_ins_overlap_parents[ 'is_close' ] == True) ][ 3 ].to_list()
    else:
        ins_names_overlap_parents = [ '' ]

    print('INS in parents: ', ins_names_overlap_parents)
    #print("Number of insertions supported by raw evidence in parents:", str(len(ins_names_overlap_parents)))

    #ins_names_overlap = ins_names_overlap_proband not in ins_names_overlap_parents
    ins_names_overlap = [x for x in ins_names_overlap_proband if x not in ins_names_overlap_parents]
    print(ins_names_overlap)
    #print(ins_names_overlap)
    exit()











    if len(ins_names_overlap) != 0:
        bed_filt_ins_new = bed_filt_ins_ref.filter(lambda x: x.name in ins_names_overlap).saveas() # create bedtool object of all SVs in ins_names_overlap
    #print(bed_filt_ins_new)
    #print(bed_filt_ins_new[1].fields[5])
    
        print(ins_names_overlap)
        for i in range(len(bed_filt_ins_new)): #for each SV
            parent_list = getParents(bed_filt_ins_new[i].fields[5], ped) #get parents
            start = int(bed_filt_ins_new[i].fields[2]) - 500 #get start and end position +- 100
            end = int(bed_filt_ins_new[i].fields[2]) + 500
            check_parents = raw_bed_ref.filter(lambda b: b.score in parent_list and b.start > start and b.start < end).saveas().to_dataframe(header=None) #check if there is a parent in raw evidence that has SV with the same start and end
        #print(check_parents)
        #if check_parents[1].fields[5] == parent_list[1] | parent_list[2]:
            #ins_names_overlap.remove(bed_filt_ins_ref[i].fields[5])
            if (bed_filt_ins_new[i].fields[3] not in ins_names_overlap):
                ins_names_overlap.append(bed_filt_ins_new[i].fields[3]) # if it is not in the list that means that it has been removed before becasuse there were multiple with this name, so readd it and check if it needs to be deleted
            if check_parents.shape[0] != 0: #if there was a parent with similar SV
            #print(bed_filt_ins_new[i].fields[3])
            #print(ins_names_overlap)
                ins_names_overlap.remove(bed_filt_ins_new[i].fields[3]) #remove SV from ins_names_overlap because it is inherited
        print(ins_names_overlap)
        exit()

       

    ## Large CNVS: Reciprocal overlap >0.4%
    verbosePrint('Checking large cnvs in raw files', verbose)
    large_bed_filt_cnv = bed_filt[bed_filt['SVTYPE'].isin(['DEL', 'DUP']) & bed_filt['is_large_cnv'] == True]

    print("Number of large CNVs found:",str(len(large_bed_filt_cnv)))

    if (len(keep_large) != 0):
        if (len(large_bed_filt_cnv.index) > 0):
            bed_filt_cnv_ref = large_bed_filt_cnv[cols_keep2].to_string(header=False, index=False)
            bed_filt_cnv_ref = pybedtools.BedTool(bed_filt_cnv_ref, from_string=True).sort()
            large_bed_filt_cnv_overlap = bed_filt_cnv_ref.intersect(raw_bed_ref,
                                                            wo=True,
                                                            f=large_raw_overlap,
                                                            r=True
                                                            ).to_dataframe(disable_auto_names=True, header=None)
            if (len(large_bed_filt_cnv_overlap) != 0):
                large_cnv_names_overlap = large_bed_filt_cnv_overlap[3].to_list()
       
            else:
                large_cnv_names_overlap = ['']
    else:
        large_cnv_names_overlap = ['']

    print("Number of large CNVs supported by raw evidence:",str(len(large_cnv_names_overlap)))


    ## Small CNVs - Reciprocal overlap >0.8%
    verbosePrint('Checking small cnvs in raw files', verbose)

    small_bed_filt_cnv = bed_filt[bed_filt['SVTYPE'].isin(['DEL', 'DUP']) & bed_filt['is_large_cnv'] == False]
    print("Number of small CNVs found:",str(len(small_bed_filt_cnv)))

    if (len(keep_small) != 0):
        if (len(small_bed_filt_cnv.index) > 0):
            bed_filt_cnv_ref = small_bed_filt_cnv[cols_keep2].to_string(header=False, index=False)
            bed_filt_cnv_ref = pybedtools.BedTool(bed_filt_cnv_ref, from_string=True).sort()
            small_bed_filt_cnv_overlap = bed_filt_cnv_ref.intersect(raw_bed_ref,
                                                            wo=True,
                                                            f=small_raw_overlap,
                                                            r=True
                                                            ).to_dataframe(disable_auto_names=True, header=None)
            if (len(small_bed_filt_cnv_overlap) != 0):
                small_cnv_names_overlap = small_bed_filt_cnv_overlap[3].to_list()
            #if 'phase2_DEL_chr19_484' in cnv_names_overlap:
            #print('exists')
            else:
                small_cnv_names_overlap = ['']
    else:
        small_cnv_names_overlap = ['']
    print("Number of small CNVs supported by raw evidence:",str(len(small_cnv_names_overlap)))
    

    ##Filtering out INS and CNV with no raw evidence
    verbosePrint('Filtering out variants with no raw evidence', verbose)

    #bed_test = bed_filt[ (bed_filt['SVTYPE'].isin(['DEL', 'DUP', 'INS'])) |
                          #(bed_filt['name'].isin(ins_names_overlap + cnv_names_overlap)) ]

    bed_final = bed_filt[ (~bed_filt['SVTYPE'].isin(['DEL', 'DUP', 'INS'])) |
                          (bed_filt['name'].isin(ins_names_overlap + large_cnv_names_overlap + small_cnv_names_overlap)) ]

    
    ##Keep samples and outliers in sepparate files
    output = bed_final[(bed_final['sample'].isin(samples_keep))]


    output_outliers = bed_final[(~bed_final['sample'].isin(samples_keep))]

    # TO DO
    # Go back to BND and mCNV and chrY
    # When complex, test the segments of the calls?

    # Write output
    output.to_csv(path_or_buf=out_file, mode='a', index=False, sep='\t', header=True)
    output_outliers.to_csv(path_or_buf=outlier_samples, mode='a', index=False, sep='\t', header=True)

if __name__ == '__main__':
    main()