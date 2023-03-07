#!/usr/bin/env python

"""
Author: Alba Sanchis-Juan and Nicole Calamari
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
import collections
from subprocess import Popen, PIPE
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
    father = ped[(ped['IndividualID'] == sample)]['FatherID'].values[0]
    mother = ped[(ped['IndividualID'] == sample)]['MotherID'].values[0]
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
    fam = ped[(ped['IndividualID'] == sample)]['FamID'].values
    chr = row['chrom']
    if fam.size != 0:
        return ('_'.join([str(fam[0]), chr]))

def getInfoParent(row, ped, vcf, parent, field):    
    row_name = row['name']
    sample = row['sample']
    parent_id = ped[(ped['IndividualID'] == sample)][parent].values

    filt_pos = vcf[(vcf['ID'] == row_name)]['FORMAT'].str.split(':').tolist()

    if (field in filt_pos[0]):
        idx = filt_pos[0].index(field)
        parent_info = vcf[(vcf['ID'] == row_name)][parent_id[0]].str.split(':').tolist()
        parent_info_field = parent_info[0][idx]
    else:
        return('NA')

    return(parent_info_field)

def getParents(proband,ped):
    mother = ped[(ped['IndividualID'] == proband)]['MotherID'].values
    father = ped[(ped['IndividualID'] == proband)]['FatherID'].values
    parent_list = [mother,father]
    return(parent_list)

def getFamilyID(row,ped):
    samp = row['sample']
    family_id = ped[(ped['IndividualID'] == samp)]['FamID'].values
    return(family_id[0])

def getBincovMatrixUrl(sample,sample_batches,batch_bincov):
    batch = sample_batches.loc[sample_batches['sample'] == sample]['batch'].iloc[0]
    bincov_gs = batch_bincov.loc[batch_bincov['batch'] == batch]['bincov'].iloc[0]
    bincov = bincov_gs.replace('gs://', '/cromwell_root/')
    return(bincov)

def writeToFilterFile(file,header,original_df,filtered_df):
    file.write(header)
    file.write("\n")
    original = original_df['name_famid']
    filtered = filtered_df['name_famid']
    to_write = list(set(original) - set(filtered))
    file.write(str(to_write))
    file.write("\n")

def writeToSizeFile(file,header,df):
    file.write(header)
    file.write(str(len(df)))
    file.write("\n")

def minGQ(row, gq_min):
    paternal_srgq = row['paternal_srgq']
    maternal_srgq = row['maternal_srgq']
    paternal_pegq = row['paternal_pegq']
    maternal_pegq = row['maternal_pegq']
    paternal_gq = row['paternal_gq']
    maternal_gq = row['maternal_gq']
    minimum = min([paternal_srgq, maternal_srgq, paternal_gq, maternal_gq, paternal_pegq, maternal_pegq])
    if minimum != '.':
        if int(min([paternal_srgq, maternal_srgq, paternal_gq, maternal_gq, paternal_pegq, maternal_pegq])) <= gq_min:
            return('Remove')
        else:
            return('Keep')
    else:
        return('Keep')
    
def convertToBedtool(df,cols_to_keep=None,sort=True):
    if cols_to_keep == None:
        string_obj = df.to_string(header=False, index=False)
        bt_obj = pybedtools.BedTool(string_obj, from_string=True)
        if sort == True:
            string_obj = df.to_string(header=False, index=False)
            bt_obj = pybedtools.BedTool(string_obj, from_string=True).sort()
    else:
        string_obj = df[cols_to_keep].to_string(header=False, index=False)
        bt_obj = pybedtools.BedTool(string_obj, from_string=True).sort()
    return(bt_obj)

def tabix_query(filename, chrom, start, end):
    #Call tabix and generate an array of strings for each line it returns.
        query = '{}:{}-{}'.format(chrom, start, end)
        process = Popen(['tabix', '-h', filename, query], stdout=PIPE)
        for line in process.stdout:
            yield line.strip().split()

def getPerSampleMatrix(row,ped,sample_batches,batch_bincov,family_member):
    chrom = row['chrom']
    start = int(row['start'])
    end = int(row['end'])
    proband = row['sample']
    proband_matrix = getBincovMatrixUrl(proband,sample_batches,batch_bincov)
    mom = getParents(proband,ped)[0][0]
    mom_matrix = getBincovMatrixUrl(mom, sample_batches, batch_bincov)
    dad = getParents(proband,ped)[1][0]
    dad_matrix = getBincovMatrixUrl(dad, sample_batches, batch_bincov)
    if family_member == 'proband':
        return([proband,proband_matrix])
    if family_member == 'mother':
        return([mom,mom_matrix])
    if family_member == 'father':
        return([dad,dad_matrix])
    
def findCoverage(row,ped,sample_batches,batch_bincov,family_member):
    chrom = row['chrom']
    start = int(row['start'])
    end = int(row['end'])
    #if family_member == 'proband':
        #matrix = getPerSampleMatrix(row,ped,sample_batches,batch_bincov,family_member='proband')[1]
    if family_member == 'mother':
        matrix = getPerSampleMatrix(row,ped,sample_batches,batch_bincov,family_member='mother')[1]
    if family_member == 'father':
        matrix = getPerSampleMatrix(row,ped,sample_batches,batch_bincov,family_member='father')[1]
    chunks = tabix_query(matrix, chrom, start, end)
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
    #if family_member == 'proband':
        #indv = getPerSampleMatrix(row,ped,sample_batches,batch_bincov,family_member='proband')[0]
    if family_member == 'mother':
        indv = getPerSampleMatrix(row,ped,sample_batches,batch_bincov,family_member='mother')[0]
    if family_member == 'father':
        indv = getPerSampleMatrix(row,ped,sample_batches,batch_bincov,family_member='father')[0]
    sample_cov_matrix = cov_matrix[indv].tolist()
    sample_cov_matrix_decoded = []
    for x in sample_cov_matrix:
        x = x.decode()
        sample_cov_matrix_decoded.append(x)
    return(sample_cov_matrix_decoded)

def getMedianCoverage(mom_matrix,dad_matrix,coverage_cutoff):
    mom_coverage = mom_matrix
    dad_coverage = dad_matrix
        
    coverage_d = {
    'Mother' : mom_coverage,
    'Father': dad_coverage
    }

    coverage_df = pd.DataFrame(coverage_d)
    median = coverage_df.median()
    median_list = median.tolist()
    median_filtered = [x for x in median_list if x < coverage_cutoff]
    if len(median_filtered) > 0:
        return('Remove')
    else:
        return('Keep')

def getCnvIntersectionDepth(bed,raw,overlap):
    intersect = bed.coverage(raw).to_dataframe(disable_auto_names=True, header=None)
    if (len(intersect) != 0):
        names_overlap = intersect[intersect[10] > overlap][6].to_list()
    else:
        names_overlap = ['']
    return(names_overlap)

def getCnvIntersectionOther(bed,raw,overlap):
    overlap = bed.intersect(raw, wo=True, f=overlap, r=True).to_dataframe(disable_auto_names=True, header=None)
    if (len(overlap) != 0):
        names_overlap = overlap[6].to_list()
    else:
        names_overlap = ['']
    return(names_overlap)

def getInsertionIntersection(bed, raw):
    overlap = bed.closest(raw).to_dataframe(disable_auto_names=True, header=None)
    if (len(overlap) != 0):
        overlap['is_close'] = abs(overlap[8] - overlap[1]) < 100
        if (sum(bool(x) for x in overlap['is_close']) > 0):
            names_overlap = overlap[(overlap['is_close'] == True)][6].to_list()
        else:
            names_overlap = ['']
    return(names_overlap)
        
def main():
    """
    Parse arguments
vcf_metrics    """
    parser = argparse.ArgumentParser(description='Parse arguments')
    parser.add_argument('--bed', dest='bed', help='Input BED file')
    parser.add_argument('--ped', dest='ped', help='Ped file')
    parser.add_argument('--vcf', dest='vcf', help='VCF file')
    parser.add_argument('--disorder', dest='disorder', help='Genomic disorder regions')
    parser.add_argument('--out', dest='out', help='Output file')
    parser.add_argument('--raw_proband', dest='raw_proband', help='Directory with raw SV calls - output from m04')
    parser.add_argument('--raw_parents', dest='raw_parents', help='Directory with raw SV calls - output from m04')
    parser.add_argument('--raw_depth_proband', dest='raw_depth_proband', help='Directory with raw SV depth calls - output from m04')
    parser.add_argument('--raw_depth_parents', dest='raw_depth_parents', help='Directory with raw depth SV calls - output from m04')
    parser.add_argument('--outliers', dest='outliers', help='Output file with SV calls in outlier samples')
    parser.add_argument('--config', dest='config', help='Config file')
    parser.add_argument('--filtered', dest='filtered_out', help='Output of filtered out variants file')
    parser.add_argument('--size_file', dest='size_file', help='Output of size of remaining bed file')
    parser.add_argument('--exclude_regions', dest='exclude_regions', help='File containing regions with known somatic mutations')
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
    raw_file_depth_proband = args.raw_depth_proband
    raw_file_depth_parent = args.raw_depth_parents
    outlier_samples = args.outliers
    verbose = args.verbose
    config_file = args.config
    filtered_out = args.filtered_out
    size = args.size_file
    exclude_regions = args.exclude_regions
    coverage = args.coverage
    batches = args.sample_batches

    with open(config_file, "r") as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    large_cnv_size = int(config['large_cnv_size'])
    gnomad_col = config['gnomad_col']
    alt_gnomad_col = config['alt_gnomad_col']
    gnomad_AF = float(config['gnomad_AF'])
    parents_AF = float(config['parents_AF'])
    large_raw_overlap = float(config['large_raw_overlap'])
    small_raw_overlap = float(config['small_raw_overlap'])
    cohort_AF = float(config['cohort_AF'])
    coverage_cutoff = float(config['coverage_cutoff'])
    depth_only_size = float(config['depth_only_size'])
    parents_overlap = float(config['parents_overlap'])
    gq_min = float(config['gq_min'])

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
    raw_bed_depth_child = pd.read_csv(raw_file_depth_proband, sep='\t', names= raw_bed_colnames, header=None).replace(np.nan, '', regex=True)
    raw_bed_depth_parent = pd.read_csv(raw_file_depth_parent, sep='\t', names= raw_bed_colnames, header=None).replace(np.nan, '', regex=True)
    exclude_regions = pd.read_csv(exclude_regions, sep='\t').replace(np.nan, '', regex=True)
    bincov_colnames = colnames=['batch', 'bincov'] 
    sample_batches_colnames = colnames=['sample', 'batch'] 
    bincov = pd.read_csv(coverage, sep='\t', names= bincov_colnames, header=None).replace(np.nan, '', regex=True)
    sample_batches = pd.read_csv(batches, sep='\t', names= sample_batches_colnames, header=None).replace(np.nan, '', regex=True)

    # Open files for writing
    filtered_out_file = open(filtered_out, "w")
    size_file = open(size, "w")

    # Get parents and children ids
    verbosePrint('Getting parents/children/affected/unaffected IDs', verbose)
    list_of_trios = [item for item, count in collections.Counter(ped["FamID"]).items() if count == 3]
    families = ped[(ped['FamID'].isin(list_of_trios)) &
                   (ped['FatherID'] != "0") &
                   (ped['MotherID'] != "0")]['FamID'].values
    trios = ped[(ped['FamID'].isin(families))]
    parents = trios[(trios['FatherID'] == "0") & (trios['MotherID'] == "0")]['IndividualID'].values
    children = trios[(trios['FatherID'] != "0") & (trios['MotherID'] != "0")]['IndividualID'].values

    # Get counts at cohort level and number of parents/children
    verbosePrint('Getting counts', verbose)
    bed['num_children'] = bed.apply(lambda r: getCount(r, children), axis=1)
    bed['num_parents'] = bed.apply(lambda r: getCount(r, parents), axis=1)
    bed['AF_parents'] = bed.apply(lambda r: getParentsFrequency(r, parents), axis=1)

    # Remove mCNVs, BNDs and SVs in sex chromosomes
    verbosePrint('Remove BND and mCNV', verbose)
    bed = bed[(~bed['svtype'].isin(['BND', 'CNV']))]
    # bed = bed[(~bed['svtype'].isin(['BND', 'CNV']) & (~bed['chrom'].isin(["chrY", "chrX"])))]
    
    # Flag if small or large CNV based on large_cnv_size cutoff and flag for removal based on size for depth only calls
    verbosePrint('Flagging calls depending on size', verbose)
    bed['is_large_cnv'] = (bed['SVLEN'] >= large_cnv_size) & ((bed['svtype'] == 'DEL') | (bed['svtype'] == 'DUP'))
    bed['is_small_cnv'] = (bed['SVLEN'] < large_cnv_size) & ((bed['svtype'] == 'DEL') | (bed['svtype'] == 'DUP'))
    bed['is_depth_only'] = (bed['EVIDENCE'] == "RD")
    bed['is_depth_only_small'] = (bed['svtype'] == "DUP") & (bed['ALGORITHMS'] == "depth") & (bed['SVLEN'] <= depth_only_size)

    # Split into one row per sample
    verbosePrint('Split into one row per sample', verbose)
    bed_split = bed.assign(sample=bed.samples.str.split(",")).explode('sample')
   
    # Sepparate variants in children and parents
    verbosePrint('Sepparate variants in children and parents', verbose)
    bed_child = bed_split[bed_split['sample'].str.contains("|".join(children))]
    bed_parents = bed_split[bed_split['sample'].str.contains("|".join(parents))]
    writeToSizeFile(size_file,"Size of bed_child: ",bed_child)

    # Flag variants in genomic disorder (GD) regions
    verbosePrint('Flagging variants in GD', verbose)
    bed_child['in_gd'] = bed_child['name'].isin(disorder[0])

    # Annotate family information for filtering
    bed_child['family_id'] = bed_child.apply(lambda r: getFamilyID(r,ped), axis=1)
    bed_child['name_famid'] = bed_child['name'] + "_" + bed_child['family_id'].astype(str).str.strip("[]")

    # Filter out by frequency - AF gnomad < 0.01 OR inGD
    verbosePrint('Filtering by frequency', verbose)
    bed_child["AF"] = pd.to_numeric(bed_child["AF"])
    try:
        bed_child[gnomad_col] = pd.to_numeric(bed_child[gnomad_col])
        bed_child_tmp = bed_child[( ((bed_child[gnomad_col] <= gnomad_AF) & (bed_child['AF'] <= cohort_AF)) |
                           ((bed_child[gnomad_col].isnull() )  & (bed_child['AF'] <= cohort_AF)) |
                          (bed_child['in_gd'] == True) )]  
    except KeyError:
        bed_child[alt_gnomad_col] = pd.to_numeric(bed_child[alt_gnomad_col])
        bed_child_tmp = bed_child[( ((bed_child[alt_gnomad_col] <= gnomad_AF) & (bed_child['AF'] <= cohort_AF)) |
                           ((bed_child[alt_gnomad_col].isnull() )  & (bed_child['AF'] <= cohort_AF)) |
                          (bed_child['in_gd'] == True) )]  
    
    writeToFilterFile(filtered_out_file,"Removed after filtering by frequency: ",bed_child,bed_child_tmp)
    bed_child = bed_child_tmp
    writeToSizeFile(size_file,"Size of bed_child after filtering by frequency: ",bed_child)

    # Get counts within family and remove if SV in parents
    verbosePrint('Keep variants in children only', verbose)
    try:
        parents_SC = int(config['parents_SC'])
    except KeyError:
        parents_SC = len(parents)
    bed_child['num_parents_family'] = bed_child.apply(lambda r: getFamilyCount(r, ped), axis=1)
    bed_child_tmp = bed_child[ (bed_child['num_parents_family'] == 0) & #there are no parents of the affected child in the samples column
                           (bed_child['num_children'] >= 1) &   #for each line there is at least one child in the samples column
                           (bed_child['AF_parents'] <= parents_AF) & #for each line, the frequency of parents in the samples column must be < 0.05
                           (bed_child['num_parents'] <= parents_SC) ] #for each line, the number of parents must be <= parents_SC or the total number of parents in the ped
    writeToFilterFile(filtered_out_file,"Removed because of cohort wide parental support: ",bed_child,bed_child_tmp)
    bed_child = bed_child_tmp  
    writeToSizeFile(size_file,"Size of bed_child after removing because of cohort wide parental support: ",bed_child)
    
    # Extract info from the VCF file
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
    bed_child_tmp = bed_child[~((bed_child['ALGORITHMS'] == "wham") & (bed_child['GQ'] == '1'))]
    writeToFilterFile(filtered_out_file,"Removed wham only and GT=1 calls: ",bed_child,bed_child_tmp)
    bed_child = bed_child_tmp  
    writeToSizeFile(size_file,"Size of bed_child after removing wham only and GT=1 calls: ",bed_child)

    # Annotate parental information
    bed_child['paternal_gq'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'FatherID', 'GQ'), axis=1)
    bed_child['maternal_gq'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'MotherID', 'GQ'), axis=1)
    bed_child['paternal_rdcn'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'FatherID', 'RD_CN'), axis=1)
    bed_child['maternal_rdcn'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'MotherID', 'RD_CN'), axis=1)
    bed_child['paternal_pegq'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'FatherID', 'PE_GQ'), axis=1)
    bed_child['maternal_pegq'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'MotherID', 'PE_GQ'), axis=1)
    bed_child['paternal_srgq'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'FatherID', 'SR_GQ'), axis=1)
    bed_child['maternal_srgq'] = bed_child.apply(lambda r: getInfoParent(r, ped, vcf, 'MotherID', 'SR_GQ'), axis=1)

    # LARGE CNV: Check for false negative in parents: check depth in parents, independently of the calls
    verbosePrint('Large CNVs check', verbose)
    # 1. Strip out if same CN in parents and proband if not in chrX
    bed_child_tmp = bed_child.loc[((bed_child['SVLEN'] <= 5000) | ((bed_child['RD_CN'] != bed_child['maternal_rdcn']) & (bed_child['RD_CN'] != bed_child['paternal_rdcn']))) | (bed_child['chrom'] == 'chrX')]
    writeToFilterFile(filtered_out_file,"Removed if RD_CN field is same as parent and not in chrX: ",bed_child,bed_child_tmp)
    bed_child = bed_child_tmp  
    writeToSizeFile(size_file,"Size of bed_child after removing if RD_CN field is same as parent and not in chrX: ",bed_child)

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
    names_overlap = bed_overlap[(bed_overlap[9] >= parents_overlap)][3].to_list()

    # Add if overlap to bed child table
    bed_child['overlap_parent'] = (bed_child['name'].isin(names_overlap))

    # Small calls:
    # If RD,SR and < large_cnv_size, treat RD,SR as SR
    verbosePrint('Small CNVs check', verbose)
    bed_child['EVIDENCE_FIX'] = bed_child['EVIDENCE']
    bed_child[(bed_child['SVLEN'] <= large_cnv_size) & \
                    (bed_child['EVIDENCE'] == "RD,SR") & \
                    ((bed_child['svtype'] == 'DEL') | (bed_child['svtype'] == 'DUP'))]['EVIDENCE_FIX'] = "SR"

    # Add if variant contains RD evidence
    bed_child['contains_RD'] = bed_child.EVIDENCE_FIX.str.contains('RD')
    
    ###################
    ##RAW FILES CHECK##
    ###################
    # Keep all CPT and CTX - they are rare and the type is different than in the raw files
    # Check only DEL,DUP and INS for raw evidence
    verbosePrint('Checking raw files', verbose)
    # Remove depth raw calls > 1MB
    raw_bed_depth_parent['SVLEN'] = raw_bed_depth_parent['end'] - raw_bed_depth_parent['start']
    raw_bed_ref_depth_parent_subset = raw_bed_depth_parent[(raw_bed_depth_parent['SVLEN'] < 1000000)]

    # Reformat raw files
    raw_bed_ref_child = convertToBedtool(raw_bed_child, sort=False)
    raw_bed_ref_parent = convertToBedtool(raw_bed_parent, sort=False)
    raw_bed_ref_depth_child = convertToBedtool(raw_bed_depth_child, sort=False)
    raw_bed_ref_depth_parent = convertToBedtool(raw_bed_ref_depth_parent_subset, sort=False)

    #Reformat de novo filt calls
    bed_child['chrom_type_sample'] = bed_child['chrom'] + "_" + bed_child['SVTYPE'] + "_" + bed_child['sample']
    bed_child['chrom_type_family'] = bed_child['chrom'] + "_" + bed_child['SVTYPE'] + "_" + bed_child['family_id'].astype(str)

    cols_keep_child = ['chrom_type_sample', 'start', 'end', 'name', 'svtype', 'sample', 'name_famid']
    cols_keep_parent = ['chrom_type_family', 'start', 'end', 'name', 'svtype', 'sample', 'name_famid']

    ##INSERTIONS:
    verbosePrint('Checking insertions in raw files', verbose)
    bed_filt_ins = bed_child[bed_child['SVTYPE'] == 'INS']
    if (len(bed_filt_ins.index) > 0):
        verbosePrint('Checking if insertion in proband is in raw files', verbose)
        bed_filt_ins_proband = convertToBedtool(bed_filt_ins, cols_to_keep=cols_keep_child,sort=True)
        ins_names_overlap_proband = getInsertionIntersection(bed_filt_ins_proband, raw_bed_ref_child)
        bed_filt_ins_tmp = bed_filt_ins[bed_filt_ins['name_famid'].isin(ins_names_overlap_proband)]
        writeToFilterFile(filtered_out_file,"Removed if insertion is not supported by raw evidence: ",bed_filt_ins, bed_filt_ins_tmp)
        verbosePrint('Checking if insertion in proband are also in raw files for the parents', verbose)
        bed_filt_ins_fam = convertToBedtool(bed_filt_ins, cols_to_keep=cols_keep_parent,sort=True)
        ins_names_overlap_parent = getInsertionIntersection(bed_filt_ins_fam, raw_bed_ref_parent)
        bed_filt_ins_tmp = bed_filt_ins[~(bed_filt_ins['name_famid'].isin(ins_names_overlap_parent))]
        writeToFilterFile(filtered_out_file,"Removed if insertion is supported by raw evidence in parent: ",bed_filt_ins, bed_filt_ins_tmp)
        ins_names_overlap = [x for x in ins_names_overlap_proband if x not in ins_names_overlap_parent]
    else:
        ins_names_overlap = ['']
    
    ## Large CNVS: Reciprocal overlap large_raw_overlap
    verbosePrint('Checking large cnvs in raw files', verbose)
    large_bed_filt_cnv_depth = bed_child[(bed_child['is_large_cnv'] == True) & (bed_child['SVLEN'] >= 5000)]
    large_bed_filt_cnv_other = bed_child[(bed_child['is_large_cnv'] == True) & (bed_child['SVLEN'] < 5000)]
    if (len(large_bed_filt_cnv_other.index) > 0):
        verbosePrint('Checking if intermediate cnv in proband is in raw files', verbose)
        bed_filt_cnv_proband_other = convertToBedtool(large_bed_filt_cnv_other, cols_to_keep=cols_keep_child,sort=True)
        large_cnv_names_overlap_proband_other = getCnvIntersectionOther(bed_filt_cnv_proband_other,raw_bed_ref_child,large_raw_overlap)
        large_bed_filt_cnv_tmp_other = large_bed_filt_cnv_other[large_bed_filt_cnv_other['name_famid'].isin(large_cnv_names_overlap_proband_other)]
        writeToFilterFile(filtered_out_file,"Removed if intermediate cnv is not supported by raw evidence: ",large_bed_filt_cnv_other, large_bed_filt_cnv_tmp_other)
        verbosePrint('Checking if intermediate cnvs in proband are also in raw files for the parents', verbose)
        bed_filt_cnv_fam_other = convertToBedtool(large_bed_filt_cnv_other, cols_to_keep=cols_keep_parent,sort=True)
        large_cnv_names_overlap_parent_other = getCnvIntersectionOther(bed_filt_cnv_fam_other,raw_bed_ref_parent,large_raw_overlap)
        large_bed_filt_cnv_tmp_other = large_bed_filt_cnv_other[~(large_bed_filt_cnv_other['name_famid'].isin(large_cnv_names_overlap_parent_other))]
        writeToFilterFile(filtered_out_file,"Removed if intermediate cnv is supported by raw evidence in parent: ",large_bed_filt_cnv_other, large_bed_filt_cnv_tmp_other)
        large_cnv_names_overlap_other = [x for x in large_cnv_names_overlap_proband_other if x not in large_cnv_names_overlap_parent_other]
    else:
        large_cnv_names_overlap_other = ['']
    if (len(large_bed_filt_cnv_depth.index) > 0):
        verbosePrint('Checking if large cnv in proband is in raw files', verbose)
        bed_filt_cnv_proband_depth = convertToBedtool(large_bed_filt_cnv_depth, cols_to_keep=cols_keep_child,sort=True)
        large_cnv_names_overlap_proband_depth = getCnvIntersectionDepth(bed_filt_cnv_proband_depth,raw_bed_ref_depth_child,large_raw_overlap)
        large_bed_filt_cnv_tmp_depth = large_bed_filt_cnv_depth[large_bed_filt_cnv_depth['name_famid'].isin(large_cnv_names_overlap_proband_depth)]
        writeToFilterFile(filtered_out_file,"Removed if large cnv is not supported by raw evidence: ",large_bed_filt_cnv_depth, large_bed_filt_cnv_tmp_depth)
        verbosePrint('Checking if large cnv in proband are also in raw files for the parents', verbose)
        bed_filt_cnv_fam_depth = convertToBedtool(large_bed_filt_cnv_depth, cols_to_keep=cols_keep_parent,sort=True)
        large_cnv_names_overlap_parent_depth = getCnvIntersectionDepth(bed_filt_cnv_fam_depth,raw_bed_ref_depth_parent,large_raw_overlap)
        large_bed_filt_cnv_tmp_depth = large_bed_filt_cnv_depth[~(large_bed_filt_cnv_depth['name_famid'].isin(large_cnv_names_overlap_parent_depth))]
        writeToFilterFile(filtered_out_file,"Removed if large cnv is supported by raw evidence in parent: ",large_bed_filt_cnv_depth, large_bed_filt_cnv_tmp_depth)
        large_cnv_names_overlap_depth = [x for x in large_cnv_names_overlap_proband_depth if x not in large_cnv_names_overlap_parent_depth]
    else:
        large_cnv_names_overlap_depth = ['']
    large_cnv_names_overlap = large_cnv_names_overlap_other + large_cnv_names_overlap_depth

    ## Small CNVs - Reciprocal overlap small_raw_overlap
    verbosePrint('Checking small cnvs in raw files', verbose)
    #small_bed_filt_cnv = bed_child[bed_child['is_small_cnv'] == True]
    #small_bed_filt_cnv_depth = bed_child[(bed_child['is_small_cnv'] == True) & (bed_child['ALGORITHMS'] == "depth") & (bed['SVLEN'] >= 5000)]
    small_bed_filt_cnv = bed_child[(bed_child['is_small_cnv'] == True)] # we do not want to check against raw depth files if CNV <= 5kb
    if (len(small_bed_filt_cnv.index) > 0):
        verbosePrint('Checking if small cnv in proband is in raw files', verbose)
        bed_filt_cnv_proband = convertToBedtool(small_bed_filt_cnv, cols_to_keep=cols_keep_child,sort=True)
        small_cnv_names_overlap_proband = getCnvIntersectionOther(bed_filt_cnv_proband,raw_bed_ref_child,small_raw_overlap)
        small_bed_filt_cnv_tmp = small_bed_filt_cnv[small_bed_filt_cnv['name_famid'].isin(small_cnv_names_overlap_proband)]
        writeToFilterFile(filtered_out_file,"Removed if small cnv is not supported by raw evidence: ",small_bed_filt_cnv, small_bed_filt_cnv_tmp)
        verbosePrint('Checking small cnvs in probands are also in raw files for parents', verbose)
        bed_filt_cnv_fam = convertToBedtool(small_bed_filt_cnv, cols_to_keep=cols_keep_parent,sort=True)
        small_cnv_names_overlap_parent = getCnvIntersectionOther(bed_filt_cnv_fam,raw_bed_ref_parent,small_raw_overlap)
        small_bed_filt_cnv_tmp = small_bed_filt_cnv[~(small_bed_filt_cnv['name_famid'].isin(small_cnv_names_overlap_parent))]
        writeToFilterFile(filtered_out_file,"Removed if small cnv is supported by raw evidence in parent: ",small_bed_filt_cnv, small_bed_filt_cnv_tmp)
        small_cnv_names_overlap = [x for x in small_cnv_names_overlap_proband if x not in small_cnv_names_overlap_parent]
    else:
        small_cnv_names_overlap = ['']
    
    bed_child = bed_child[(~bed_child['SVTYPE'].isin(['DEL', 'DUP', 'INS'])) | (bed_child['name_famid'].isin(ins_names_overlap + large_cnv_names_overlap + small_cnv_names_overlap))]
    writeToSizeFile(size_file,"Size of bed_child after checking raw files: ",bed_child)


    ###############
    ## FILTERING ##
    ###############
    verbosePrint('Filtering out calls', verbose)
    # 1. Filter by region
    #  Keep if in GD region
    keep_gd = bed_child[(bed_child['in_gd'] == True)]['name_famid'].to_list()
    
    # Filter out calls in exclude regions
    # Reformat exclude_regions to bedtool
    exclue_regions_bt = convertToBedtool(exclude_regions,sort=True)
    # convert bed_final to bedtool
    cols_keep_exclude_regions = ['chrom', 'start', 'end', 'name', 'svtype', 'sample', 'name_famid']
    bed_child_bt = convertToBedtool(bed_child,cols_keep_exclude_regions, sort=True)
    exclude_regions_intersect = bed_child_bt.coverage(exclue_regions_bt).to_dataframe(disable_auto_names=True, header=None) #HB said to use bedtools coverage, -f and -F give the same SVs to be removed
    if (len(exclude_regions_intersect) != 0):
        remove_regions = exclude_regions_intersect[exclude_regions_intersect[10] > 0.5][6].to_list()
        bed_child_tmp = bed_child[~(bed_child['name_famid'].isin(remove_regions))]
        writeToFilterFile(filtered_out_file,"Removed if overlaps with exclude region: ",bed_child,bed_child_tmp)
    else:
        remove_regions = ['']
    
    # 2. Filter by type
    # Keep any SV type that is not a DEL, DUP, or INS
    keep_other_sv = bed_child[~bed_child['SVTYPE'].isin(['DEL', 'DUP', 'INS'])]['name_famid'].to_list()

    # 3. Filter on DELs, DUPs, and INS
    # Filter by size
    # Filter out if large CNVs have parents overlap
    remove_large = bed_child[(bed_child['is_large_cnv'] == True) &
                            (bed_child['overlap_parent'] == True)]['name_famid'].to_list()
    bed_child_tmp = bed_child[~(bed_child['name_famid'].isin(remove_large))]
    writeToFilterFile(filtered_out_file,"Removed if large CNV that has parents overlap: ",bed_child,bed_child_tmp)

    # Filter out if small cnvs that are SR-only don't have BOTHSIDES_SUPPORT
    remove_small = bed_child[(bed_child['is_small_cnv'] == True) &
                             (bed_child['EVIDENCE_FIX'] == 'SR') &
                             ~(bed_child.FILTER.str.contains('BOTHSIDES_SUPPORT'))]['name_famid'].to_list()
    bed_child_tmp = bed_child[~(bed_child['name_famid'].isin(remove_small))]
    writeToFilterFile(filtered_out_file,"Removed if small cnv that is SR-only and does not have BOTHSIDES_SUPPORT: ",bed_child,bed_child_tmp)

    # Filter out calls that are depth only and < depth_only_size
    remove_depth_small = bed_child[ (bed_child['is_depth_only_small'] == True) ]['name_famid'].to_list()
    bed_child_tmp = bed_child[~(bed_child['name_famid'].isin(remove_depth_small))]
    writeToFilterFile(filtered_out_file,"Removed if depth only and < 10kb \n",bed_child,bed_child_tmp)

    # Filter out DELs that are >500bp and RD_CN=2 and PE only envidence
    remove_dels = bed_child[(bed_child['SVTYPE'] == 'DEL') & ((bed_child['RD_CN'] == '2') | (bed_child['RD_CN'] == '3')) & (bed_child['EVIDENCE'] == 'PE')]['name_famid'].to_list()
    bed_child_tmp = bed_child[~(bed_child['name_famid'].isin(remove_dels))]
    writeToFilterFile(filtered_out_file,"Removed if DEL with RD_CN >= 2 and PE only evidence: ",bed_child,bed_child_tmp)

    # 4. Filter by quality
    # Filter out if parents GQ is <= gq_min
    bed_child['keep_gq'] = bed_child.apply(lambda r: minGQ(r, gq_min), axis=1)
    remove_gq = bed_child[bed_child['keep_gq'] == 'Remove']['name_famid'].to_list()
    bed_child_tmp = bed_child[~(bed_child['name_famid'].isin(remove_gq))]
    writeToFilterFile(filtered_out_file,"Removed if minimum parental GQ is <= gq_min: ",bed_child,bed_child_tmp)

    # Remove if low coverage in parents
    bed_child['median_coverage'] = bed_child.apply(lambda r: getMedianCoverage(findCoverage(r, ped, sample_batches, bincov, family_member='mother'),findCoverage(r, ped, sample_batches, bincov, family_member='father'), coverage_cutoff), axis=1)
    remove_coverage = bed_child[bed_child['median_coverage'] == 'Remove']['name_famid'].to_list()
    bed_child_tmp = bed_child[~(bed_child['name_famid'].isin(remove_coverage))]
    writeToFilterFile(filtered_out_file,"Removed if median coverage in mother or father is <= coverage_cutoff: ",bed_child,bed_child_tmp)
    
    # 5. Clean up and remove duplicated CPX SV
    # Remove duplicated CPX events that come from vcf output of module 18
    bed_child['is_cpx'] = (bed_child['SVTYPE'] == "CPX")
    bed_child['is_duplicated'] = bed_child.duplicated(keep='first', subset=['start', 'end', 'sample'])
    remove_duplicated = bed_child[(bed_child['is_duplicated'] == True) & (bed_child['is_cpx'] == True)]['name_famid'].to_list()
    bed_child_tmp = bed_child[~(bed_child['name_famid'].isin(remove_duplicated))]
    writeToFilterFile(filtered_out_file,"Removed CPX if duplicated SV: ",bed_child,bed_child_tmp)
    
    # Join lists to remove and keep unique values
    keep_names = list(set(keep_gd + keep_other_sv))
    remove_names = list(set(remove_regions + remove_large + remove_small + remove_depth_small + remove_dels + remove_gq + remove_coverage))
    final_remove_names = [x for x in remove_names if x not in keep_names]
    
    # Subset table
    bed_filt = bed_child[~(bed_child['name_famid'].isin(final_remove_names))]
    bed_final = bed_filt[~(bed_filt['name_famid'].isin(remove_duplicated))] #remove duplicated CPX last
    writeToSizeFile(size_file,"Size of bed_filt after filtering: ",bed_final)

    #####################
    ## REMOVE OUTLIERS ##
    #####################
    verbosePrint('Removing outliers', verbose)
    sample_counts = bed_final['sample'].value_counts()
    q3, q1 = np.percentile(sample_counts.tolist(), [75, 25])
    iqr = q3 - q1
    count_threshold = int(math.ceil((1.5*iqr + q3)*3))
    sample_counts2 = sample_counts.rename_axis('sample').reset_index(name='count')
    samples_keep = sample_counts2[(sample_counts2['count'] <= count_threshold)]['sample'].tolist()

    #move sample column to column 6
    col = bed_final.pop("sample")
    bed_final.insert(loc=5, column='sample', value=col, allow_duplicates=True)

    ##Keep samples and outliers in sepparate files
    output = bed_final[(bed_final['sample'].isin(samples_keep))]
    output_outliers = bed_final[(~bed_final['sample'].isin(samples_keep))]

    # Close out files
    filtered_out_file.close()
    size_file.close()

    # TO DO
    # Go back to BND and mCNV and chrY
    # When complex, test the segments of the calls?

    # Write output
    output.to_csv(path_or_buf=out_file, mode='a', index=False, sep='\t', header=True)
    output_outliers.to_csv(path_or_buf=outlier_samples, mode='a', index=False, sep='\t', header=True)

if __name__ == '__main__':
    main()