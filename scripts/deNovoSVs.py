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
import tabix
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
    return(family_id[0])

def getBincovMatrix(sample,sample_batches,batch_bincov):
    batch = sample_batches.loc[sample_batches['sample'] == sample]['batch'].iloc[0]
    bincov_gs = batch_bincov.loc[batch_bincov['batch'] == batch]['bincov'].iloc[0]
    bincov = bincov_gs.replace('gs://', '/cromwell_root/')
    return(bincov)


    

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
    parser.add_argument('--raw_proband', dest='raw_proband', help='Directory with raw SV calls - output from m04')
    parser.add_argument('--raw_parents', dest='raw_parents', help='Directory with raw SV calls - output from m04')
    parser.add_argument('--outliers', dest='outliers', help='Output file with SV calls in outlier samples')
    parser.add_argument('--config', dest='config', help='Config file')
    parser.add_argument('--filtered', dest='filtered_out', help='Output file')
    parser.add_argument('--SM_regions', dest='somatic_mutation_regions', help='File containing regions with known somatic mutations')
    parser.add_argument('--coverage', dest='coverage', help='File with batch in first column respective coverage file in second column')
    parser.add_argument('--sample_batches', dest='sample_batches', help='File with samples in first column and their respective batch in second column')
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
    filtered_out = args.filtered_out
    somatic_mutation_regions = args.somatic_mutation_regions
    coverage = args.coverage
    batches = args.sample_batches


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
    coverage_cutoff = float(config['coverage_cutoff'])

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
    b = pd.read_csv(somatic_mutation_regions, sep='\t').replace(np.nan, '', regex=True)
    bincov_colnames = colnames=['batch', 'bincov'] 
    sample_batches_colnames = colnames=['sample', 'batch'] 
    bincov = pd.read_csv(coverage, sep='\t', names= bincov_colnames, header=None).replace(np.nan, '', regex=True)
    sample_batches = pd.read_csv(batches, sep='\t', names= sample_batches_colnames, header=None).replace(np.nan, '', regex=True)



    # Get parents and children ids
    verbosePrint('Getting parents/children/affected/unaffected IDs', verbose)
    
    families = ped[(ped['family_size'] >= 3) &
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
    bed['is_depth_only2'] = (bed['ALGORITHMS'] == "depth")


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

    bed_child['family_id'] = bed_child.apply(lambda r: getFamilyID(r,ped), axis=1)
    bed_child['name_famid'] = bed_child['name'] + "_" + bed_child['family_id'].astype(str).str.strip("[]")

    # Filter out by frequency - AF gnomad < 0.01 OR inGD
    verbosePrint('Filtering by frequency', verbose)
    bed_child[gnomad_col] = pd.to_numeric(bed_child[gnomad_col])
    bed_child["AF"] = pd.to_numeric(bed_child["AF"])


    f = open(filtered_out, "w")
    f.write("Removed after filtering by frequency \n")
    f.write(str(bed_child[~( ((bed_child[gnomad_col] <= gnomad_AF) & (bed_child['AF'] <= cohort_AF)) |
                           ((bed_child[gnomad_col].isnull() )  & (bed_child['AF'] <= cohort_AF)) |
                          (bed_child['in_gd'] == True) )]['name_famid'].tolist()))
    f.write("\n")


    bed_child = bed_child[( ((bed_child[gnomad_col] <= gnomad_AF) & (bed_child['AF'] <= cohort_AF)) |
                           ((bed_child[gnomad_col].isnull() )  & (bed_child['AF'] <= cohort_AF)) |
                          (bed_child['in_gd'] == True) )]  


    

    print("Size of bed_child after filtering by frequency:",str(len(bed_child)))
    

    # Get counts within family and remove if SV in parents
    verbosePrint('Keep variants in children only', verbose)
    bed_child['num_parents_family'] = bed_child.apply(lambda r: getFamilyCount(r, ped), axis=1)

    
    f.write("Removed after keeping variants in children only \n")
    f.write(str(bed_child[~( (bed_child['num_parents_family'] == 0) &
                           (bed_child['num_children'] >= 1) &
                           (bed_child['AF_parents'] <= parents_AF) &
                           (bed_child['num_parents'] <= parents_SC) )]['name_famid'].tolist()))
    f.write("\n")



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
    f.write("Removed wham only and GT=1 calls \n")
    f.write(str(bed_child[((bed_child['ALGORITHMS'] == "wham") & (bed_child['GQ'] == '1'))]['name_famid'].tolist()))
    f.write("\n")
    #print(bed_child[((bed_child['ALGORITHMS'] == "wham") & (bed_child['GQ'] == '1'))])
    bed_child = bed_child[~((bed_child['ALGORITHMS'] == "wham") & (bed_child['GQ'] == '1'))]

    print("Size of bed_child after removing wham only and GT=1 calls:",str(len(bed_child)))


    # LARGE CNV: Check for false negative in parents: check depth in parents, independently of the calls
    verbosePrint('Large CNVs check', verbose)
    # 1. Add parental RD_CN field and strip out if same as in proband
    bed_child['paternal_rdcn'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'paternal_id', 'RD_CN'), axis=1)
    bed_child['maternal_rdcn'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'maternal_id', 'RD_CN'), axis=1)
    bed_child['paternal_srgq'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'paternal_id', 'SR_GQ'), axis=1)
    bed_child['maternal_srgq'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'maternal_id', 'SR_GQ'), axis=1)
    bed_child['paternal_gq'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'paternal_id', 'GQ'), axis=1)
    bed_child['maternal_gq'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'maternal_id', 'GQ'), axis=1)
    bed_child['paternal_pegq'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'paternal_id', 'PE_GQ'), axis=1)
    bed_child['maternal_pegq'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'maternal_id', 'PE_GQ'), axis=1)

    bed_child = bed_child.loc[(bed_child['is_large_cnv'] == False) | ((bed_child['RD_CN'] != bed_child['maternal_rdcn']) & (bed_child['RD_CN'] != bed_child['paternal_rdcn']))]
    
    print("Size of bed_child after removing if RD_CN field is same as parent:",str(len(bed_child)))
    f.write("Removed if RD_CN field is same as parent\n")
    f.write(str(bed_child.loc[(bed_child['is_large_cnv'] == True) & ((bed_child['RD_CN'] == bed_child['maternal_rdcn']) | (bed_child['RD_CN'] == bed_child['paternal_rdcn']))]['name_famid'].tolist()))
    f.write("\n")
    

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
    #print(len(keep_gd))
    # Filter out if large CNVs don't have RD support and parents overlap
    keep_large = bed_child[ (bed_child['is_large_cnv'] == True) &
                            (bed_child['contains_RD'] == True) &
                            (bed_child['overlap_parent'] == False) ]['name'].to_list()
    #print("large_cnvs")
    #print(len(keep_large))
    # Filter out if small cnvs that are SR-only don't have BOTHSIDES_SUPPORT
    keep_small = bed_child[(bed_child['is_small_cnv'] == True) &
                           ( (bed_child['EVIDENCE_FIX'] != 'SR') |
                           ( (bed_child['EVIDENCE_FIX'] == 'SR') &
                             (bed_child.FILTER.str.contains('BOTHSIDES_SUPPORT'))))
                          ]['name'].to_list()
    #print("small_cnvs before filtering:")
    #print(bed_child[(bed_child['is_small_cnv'] == True)]["FILTER"])
    #print(bed_child[(bed_child['is_small_cnv'] == True)])
    #print("small cnvs after filtering")
    #print(keep_small)

    # Keep any other SV type
    keep_other_sv = bed_child[~bed_child['SVTYPE'].isin(['DEL', 'DUP'])]['name'].to_list()

    # Join lists and keep unique values
    keep_names = list(set(keep_gd + keep_large + keep_small + keep_other_sv))

    #Subset table
    bed_filt = bed_child[(bed_child['name'].isin(keep_names))]

    print("Size of bed_filt after removing large CNVs that dont have RD support and parents overlap and small cnvs that are SR only and dont have support on both sides:",str(len(bed_filt)))
    f.write("Removed if large CNV and does not have RD support and parent overlap and small cnvs that are SR only and dont have support of both sides\n")
    f.write(str(bed_child[(~bed_child['name'].isin(keep_names))]['name_famid'].tolist()))
    f.write("\n")
    

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

    #bed_filt['name_sample'] = bed_filt['name'] + "_" + bed_filt['sample']
    #bed_filt['family_id'] = bed_filt.apply(lambda r: getFamilyID(r,ped), axis=1)
    #bed_filt['name_famid'] = bed_filt['name'] + "_" + bed_filt['family_id'].astype(str).str.strip("[]")
    bed_filt['chrom_type_sample'] = bed_filt['chrom'] + "_" + bed_filt['SVTYPE'] + "_" + bed_filt['sample']
    bed_filt['chrom_type_family'] = bed_filt['chrom'] + "_" + bed_filt['SVTYPE'] + "_" + bed_filt['family_id'].astype(str)

    


    cols_keep2 = ['chrom_type_sample', 'start', 'end', 'name', 'svtype', 'sample', 'name_famid']
    cols_keep3 = ['chrom_type_family', 'start', 'end', 'name', 'svtype', 'sample', 'name_famid']


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
        bed_filt_ins_overlap_proband['is_close'] = abs(bed_filt_ins_overlap_proband[8] - bed_filt_ins_overlap_proband[1]) < 100

        ins_names_overlap_proband = bed_filt_ins_overlap_proband[(bed_filt_ins_overlap_proband['is_close'] == True)][6].to_list()
    else:
        ins_names_overlap_proband = ['']
    print("Number of insertions supported by raw evidence:",str(len(ins_names_overlap_proband)))
    
    f.write("Insertions removed because not supported by raw evidence \n")
    f.write(str([x for x in bed_filt_ins['name_famid'] if x not in ins_names_overlap_proband]))
    f.write("\n")
    

    #print(ins_names_overlap_proband)
    #exit()

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
        bed_filt_ins_overlap_parents[ 'is_close' ] = abs(bed_filt_ins_overlap_parents[ 8 ] - bed_filt_ins_overlap_parents[ 1 ]) < 100
        

        ins_names_overlap_parents = bed_filt_ins_overlap_parents[ (bed_filt_ins_overlap_parents[ 'is_close' ] == True) ][6].to_list()
    else:
        ins_names_overlap_parents = [ '' ]

    print('INS in parents: ', ins_names_overlap_parents)
    print("Number of insertions supported by raw evidence in parents:", str(len(ins_names_overlap_parents)))


    ins_names_overlap = [x for x in ins_names_overlap_proband if x not in ins_names_overlap_parents]
    print(ins_names_overlap)
    print("Final number of insertions in de novo output:", str(len(ins_names_overlap_parents)))
    #print(ins_names_overlap)

    f.write("Insertions removed because supported by raw evidence of parents \n")
    f.write(str(ins_names_overlap_parents))
    f.write("\n")


    ## Large CNVS: Reciprocal overlap >0.4%
    verbosePrint('Checking large cnvs in raw files', verbose)
    large_bed_filt_cnv = bed_filt[bed_filt['SVTYPE'].isin(['DEL', 'DUP']) & bed_filt['is_large_cnv'] == True]

    print("Number of large CNVs found:",str(len(large_bed_filt_cnv)))

    if (len(keep_large) != 0):
        if (len(large_bed_filt_cnv.index) > 0):
            bed_filt_cnv_proband = large_bed_filt_cnv[cols_keep2].to_string(header=False, index=False)
            bed_filt_cnv_proband = pybedtools.BedTool(bed_filt_cnv_proband, from_string=True).sort()
            large_bed_filt_cnv_overlap_proband = bed_filt_cnv_proband.intersect(raw_bed_ref_child,
                                                            wo=True,
                                                            f=large_raw_overlap,
                                                            r=True
                                                            ).to_dataframe(disable_auto_names=True, header=None)
            if (len(large_bed_filt_cnv_overlap_proband) != 0):
                large_cnv_names_overlap_proband = large_bed_filt_cnv_overlap_proband[6].to_list()
       
            else:
                large_cnv_names_overlap_proband = ['']
    else:
        large_cnv_names_overlap_proband = ['']

    print("Number of large CNVs supported by raw evidence:",str(len(large_cnv_names_overlap_proband)))
    print(large_cnv_names_overlap_proband)
    #exit()

    f.write("Large CNVs removed because not supported by raw evidence \n")
    f.write(str([x for x in large_bed_filt_cnv['name_famid'] if x not in large_cnv_names_overlap_proband]))
    f.write("\n")



    
    verbosePrint('Checking if large cnvs in proband are also in raw files for the parents', verbose)
    if (len(keep_large) != 0):
        if (len(large_bed_filt_cnv.index) > 0):
            bed_filt_cnv_fam = large_bed_filt_cnv[cols_keep3].to_string(header=False, index=False)
            bed_filt_cnv_fam = pybedtools.BedTool(bed_filt_cnv_fam, from_string=True).sort()
            large_bed_filt_cnv_overlap_parents = bed_filt_cnv_fam.intersect(raw_bed_ref_parent,
                                                            wo=True,
                                                            f=large_raw_overlap,
                                                            r=True
                                                            ).to_dataframe(disable_auto_names=True, header=None)
            if (len(large_bed_filt_cnv_overlap_parents) != 0):
                large_cnv_names_overlap_parents = large_bed_filt_cnv_overlap_parents[6].to_list()
       
            else:
                large_cnv_names_overlap_parents = ['']
    else:
        large_cnv_names_overlap_parents = ['']

    print("Number of large CNVs supported by raw evidence in parents:",str(len(large_cnv_names_overlap_parents)))
    print(large_cnv_names_overlap_parents)

    f.write("Large CNVs removed because supported by raw evidence of parents \n")
    f.write(str(large_cnv_names_overlap_parents))
    f.write("\n")
    
    

    large_cnv_names_overlap = [x for x in large_cnv_names_overlap_proband if x not in large_cnv_names_overlap_parents]
    print(large_cnv_names_overlap)
    print("Final number of large CNVs in de novo output:",str(len(large_cnv_names_overlap)))
    #exit()



    ## Small CNVs - Reciprocal overlap >0.8%
    verbosePrint('Checking small cnvs in raw files', verbose)

    small_bed_filt_cnv = bed_filt[bed_filt['SVTYPE'].isin(['DEL', 'DUP']) & bed_filt['is_large_cnv'] == False]
    print("Number of small CNVs found:",str(len(small_bed_filt_cnv)))

    if (len(keep_small) != 0):
        if (len(small_bed_filt_cnv.index) > 0):
            bed_filt_cnv_probands = small_bed_filt_cnv[cols_keep2].to_string(header=False, index=False)
            bed_filt_cnv_probands = pybedtools.BedTool(bed_filt_cnv_probands, from_string=True).sort()
            small_bed_filt_cnv_overlap_probands = bed_filt_cnv_probands.intersect(raw_bed_ref_child,
                                                            wo=True,
                                                            f=small_raw_overlap,
                                                            r=True
                                                            ).to_dataframe(disable_auto_names=True, header=None)
            if (len(small_bed_filt_cnv_overlap_probands) != 0):
                small_cnv_names_overlap_probands = small_bed_filt_cnv_overlap_probands[6].to_list()
            #if 'phase2_DEL_chr19_484' in cnv_names_overlap:
            #print('exists')
            else:
                small_cnv_names_overlap_probands = ['']
    else:
        small_cnv_names_overlap_probands = ['']
    print("Number of small CNVs supported by raw evidence:",str(len(small_cnv_names_overlap_probands)))
    print(small_cnv_names_overlap_probands)
   

    f.write("Small CNVs removed because not supported by raw evidence \n")
    f.write(str([x for x in small_bed_filt_cnv['name_famid'] if x not in small_cnv_names_overlap_probands]))
    f.write("\n")

    verbosePrint('Checking small cnvs in probands are also in raw files for parents', verbose)

    small_bed_filt_cnv = bed_filt[bed_filt['SVTYPE'].isin(['DEL', 'DUP']) & bed_filt['is_large_cnv'] == False]
    print("Number of small CNVs found:",str(len(small_bed_filt_cnv)))

    if (len(keep_small) != 0):
        if (len(small_bed_filt_cnv.index) > 0):
            bed_filt_cnv_fam = small_bed_filt_cnv[cols_keep3].to_string(header=False, index=False)
            bed_filt_cnv_fam = pybedtools.BedTool(bed_filt_cnv_fam, from_string=True).sort()
            small_bed_filt_cnv_overlap_parents = bed_filt_cnv_fam.intersect(raw_bed_ref_parent,
                                                            wo=True,
                                                            f=small_raw_overlap,
                                                            r=True
                                                            ).to_dataframe(disable_auto_names=True, header=None)
            if (len(small_bed_filt_cnv_overlap_parents) != 0):
                small_cnv_names_overlap_parents = small_bed_filt_cnv_overlap_parents[6].to_list()
            #if 'phase2_DEL_chr19_484' in cnv_names_overlap:
            #print('exists')
            else:
                small_cnv_names_overlap_parents = ['']
    else:
        small_cnv_names_overlap_parents = ['']
    print("Number of small CNVs supported by raw evidence in parents:",str(len(small_cnv_names_overlap_parents)))
    print(small_cnv_names_overlap_parents)

    f.write("Small CNVs removed because supported by raw evidence of parents \n")
    f.write(str(small_cnv_names_overlap_parents))
    f.write("\n")
    


    small_cnv_names_overlap = [x for x in small_cnv_names_overlap_probands if x not in small_cnv_names_overlap_parents]
    print(small_cnv_names_overlap)
    print("Final number of small CNVs in de novo output:",str(len(small_cnv_names_overlap)))

    

    ##Filtering out INS and CNV with no raw evidence
    verbosePrint('Filtering out variants with no raw evidence', verbose)

    #bed_test = bed_filt[ (bed_filt['SVTYPE'].isin(['DEL', 'DUP', 'INS'])) |
                          #(bed_filt['name'].isin(ins_names_overlap + cnv_names_overlap)) ]

    bed_final = bed_filt[ (~bed_filt['SVTYPE'].isin(['DEL', 'DUP', 'INS'])) |
                          (bed_filt['name_famid'].isin(ins_names_overlap + large_cnv_names_overlap + small_cnv_names_overlap)) ]

    
    ###################
    ##FINAL FILTERING##
    ###################
    

    #Remove if parents GQ is 0
    remove_gq_list = []
    minimum = ['0']
    for name_famid in bed_final['name_famid']:
        paternal_srgq = bed_final[bed_final['name_famid'] == name_famid]['paternal_srgq'].tolist()[0]
        maternal_srgq = bed_final[bed_final['name_famid'] == name_famid]['maternal_srgq'].tolist()[0]
        paternal_pegq = bed_final[bed_final['name_famid'] == name_famid]['paternal_pegq'].tolist()[0]
        maternal_pegq = bed_final[bed_final['name_famid'] == name_famid]['maternal_pegq'].tolist()[0]
        paternal_gq = bed_final[bed_final['name_famid'] == name_famid]['paternal_gq'].tolist()[0]
        maternal_gq = bed_final[bed_final['name_famid'] == name_famid]['paternal_gq'].tolist()[0]

  
        if min([paternal_srgq, maternal_srgq, paternal_gq, maternal_gq, paternal_pegq, maternal_pegq]) in minimum:
            print(name_famid)
            remove_gq_list.append(name_famid)


    bed_final = bed_final[(~bed_final['SVTYPE'].isin(['DEL', 'DUP', 'INS'])) | ~(bed_final['name_famid'].isin(remove_gq_list))]


    f.write("Removed if paternal GQ or maternal GQ is <5 \n")
    f.write(str(remove_gq_list))
    f.write("\n")
    




    #remove calls that are depth only and < 10kb
    f.write("Removed if depth only and < 10kb \n")
    f.write(str((bed_final[(bed_final['is_depth_only2'] == True) & (bed_final['SVLEN'] <= 10000)]['name_famid'].to_list())))
    f.write("\n")
    #print(bed_final[(bed_final['is_depth_only'] == True) & (bed_final['SVLEN'] < 10000)]['SVLEN'])
    #exit()
    bed_final = bed_final[(~bed_final['SVTYPE'].isin(['DEL', 'DUP', 'INS'])) | ~((bed_final['is_depth_only2'] == True) & (bed_final['SVLEN'] <= 10000))]
    #bed_final = bed_final[ ~((bed_final['is_depth_only'] == True) & (bed_final['SVLEN'] < 10000))]
    #print(bed_final)


   
    
    #remove if in somatic mutation region
    b_string = b.to_string(header=False, index=False)
    b_bt = pybedtools.BedTool(b_string, from_string=True).sort()
    

    #convert bed_final to bedtool

    cols_keep4 = ['chrom', 'start', 'end', 'name', 'svtype', 'sample', 'name_famid']

    bed_final_string = bed_final[cols_keep4].to_string(header=False, index=False)
    bed_final_bt = pybedtools.BedTool(bed_final_string, from_string=True).sort()
    #intersect_output = bed_final_bt.intersect(b_bt, wo=True, F=0.5).to_dataframe(disable_auto_names=True, header=None)
    intersect_output = bed_final_bt.coverage(b_bt).to_dataframe(disable_auto_names=True, header=None) #HB said to use bedtools coverage, -f and -F give the same SVs to be removed
    #print(intersect_output[intersect_output[3] == 'batch1.chr14.final_cleanup_CPX_chr14_32'])

    if (len(intersect_output) != 0):
        somatic_mutation_regions_overlap = intersect_output[intersect_output[10] > 0.5][6].to_list()
        f.write("Removed if overlaps with known somatic mutation region \n")
        f.write(str(somatic_mutation_regions_overlap))
        f.write("\n")
    
    
    bed_final = bed_final[(~bed_final['SVTYPE'].isin(['DEL', 'DUP', 'INS'])) | ~(bed_final['name_famid'].isin(somatic_mutation_regions_overlap))]
    #print(bed_final)




    #remove low coverage  
    from subprocess import Popen, PIPE
    def tabix_query(filename, chrom, start, end):
    #Call tabix and generate an array of strings for each line it returns.
        query = '{}:{}-{}'.format(chrom, start, end)
        process = Popen(['tabix', '-h', filename, query], stdout=PIPE)
        for line in process.stdout:
            yield line.strip().split()

    
    remove_for_coverage = []
    for name_famid in bed_final['name_famid']:
        chrom = bed_final[bed_final['name_famid'] == name_famid]['chrom'].tolist()[0]
        start = int(bed_final[bed_final['name_famid'] == name_famid]['start'].tolist()[0])
        end = int(bed_final[bed_final['name_famid'] == name_famid]['end'].tolist()[0])
        proband = bed_final[bed_final['name_famid'] == name_famid]['sample'].tolist()[0]
        proband_matrix = getBincovMatrix(proband, sample_batches, bincov)
        mom = getParents(proband,ped)[0][0]
        mom_matrix = getBincovMatrix(mom, sample_batches, bincov)
        dad = getParents(proband,ped)[1][0]
        dad_matrix = getBincovMatrix(dad, sample_batches, bincov)


        #proband
        chunks = tabix_query(proband_matrix, chrom, start, end)
        cov = []
        for chunk in chunks:
            cov.append(chunk)
        cov_matrix = pd.DataFrame(cov) #coverage matrix of only that region
        header = cov_matrix.iloc[0].to_list() #grab the first row for the header
        new_header = []
        for x in header:
            x = x.decode()
            new_header.append(x)
        cov_matrix = cov_matrix[1:] #take the data less the header row
        cov_matrix.columns = new_header
        proband_list = cov_matrix[proband].tolist()
        proband_list_decoded = []
        for x in proband_list:
            x = x.decode()
            proband_list_decoded.append(x)

        #mom
        chunks = tabix_query(mom_matrix, chrom, start, end)
        cov = []
        for chunk in chunks:
            cov.append(chunk)
        cov_matrix = pd.DataFrame(cov) #coverage matrix of only that region
        header = cov_matrix.iloc[0].to_list() #grab the first row for the header
        new_header = []
        for x in header:
            x = x.decode()
            new_header.append(x)
        cov_matrix = cov_matrix[1:] #take the data less the header row
        cov_matrix.columns = new_header
        mother_list = cov_matrix[mom].tolist()
        mother_list_decoded = []
        for x in mother_list:
            x = x.decode()
            mother_list_decoded.append(x)

        #dad
        chunks = tabix_query(dad_matrix, chrom, start, end)
        cov = []
        for chunk in chunks:
            cov.append(chunk)
        cov_matrix = pd.DataFrame(cov) #coverage matrix of only that region
        header = cov_matrix.iloc[0].to_list() #grab the first row for the header
        new_header = []
        for x in header:
            x = x.decode()
            new_header.append(x)
        cov_matrix = cov_matrix[1:] #take the data less the header row
        cov_matrix.columns = new_header
        father_list = cov_matrix[dad].tolist()
        father_list_decoded = []
        for x in father_list:
            x = x.decode()
            father_list_decoded.append(x)


        coverage_d = {
        'Proband': proband_list_decoded,
        'Mother' : mother_list_decoded,
        'Father': father_list_decoded
        }

        coverage_df = pd.DataFrame(coverage_d)

        median = coverage_df.median()


        median_filtered = [x for x in median if x < coverage_cutoff]
        if len(median_filtered) > 0:
            remove_for_coverage.append(name_famid)


    f.write("Removed if coverage is low for proband, mom, or dad \n")
    f.write(str(remove_for_coverage))
    f.write("\n")
    f.close()

    bed_final = bed_final[(~bed_final['SVTYPE'].isin(['DEL', 'DUP', 'INS'])) | ~(bed_final['name_famid'].isin(remove_for_coverage))]
    #print(bed_final)
    #exit()


    #Remove duplicated CPX events that come from vcf output of module 18
    bed_final['is_cpx'] = (bed_final['SVTYPE'] == "CPX")
    #cpx = bed_final[bed_final['is_cpx'] == True]
    print(bed_final)
    bed_final = bed_final.drop_duplicates(subset=['start', 'end', 'sample'])
    print(bed_final)



    '''
    tb = tabix.open(coverage)
    #tabix(tb, index)
    for name_famid in bed_final['name_famid']:
        chrom = str(bed_final[bed_final['name_famid'] == name_famid]['chrom'])
        start = int(bed_final[bed_final['name_famid'] == name_famid]['start'])
        end = int(bed_final[bed_final['name_famid'] == name_famid]['end'])
        tb.query(chrom, start, end)
        exit()
    

    bed_final_bedtools = bed_final[cols_keep2].to_string(header=False, index=False)
    bed_final_bedtools = pybedtools.BedTool(bed_final_bedtools, from_string=True).sort()
    c = bed_final_bedtools.coverage(bed_final_bedtools)
    print(c)
    exit()
    

    batch_bincov = pd.read_csv("tabix_test_2.bed", sep='\t', header=0).replace(np.nan, '', regex=True)
    print(batch_bincov)
    batch_bincov['median'] = batch_bincov.iloc[:, 3:].median(axis=1).to_list()
    print(batch_bincov[batch_bincov['median'] < 5]) #find the median of all rows and if its less than 5 remove any SV in that region from de novo list, but the problem is that the sample itself might not be low coverage at that region
    batch_bincov.to_csv(path_or_buf="median.txt", mode='a', index=False, sep='\t', header=True)


    print(batch_bincov[batch_bincov['__1883001__7ff04c'] < 10]['Start'])
    exit()

    low_coverage = pd.read_csv("tabix_test_after.bed", sep='\t', header=None).replace(np.nan, '', regex=True)
    low_coverage_IDs = low_coverage[118]
    print(low_coverage_IDs)
    bed_final = bed_final[ ~((bed_final['name_famid'].isin(low_coverage_IDs)))]

    print(bed_final)
    exit()
    '''

    

    
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