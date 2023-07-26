import sys,os,argparse
import numpy as np
import pandas as pd
#[_,varfile,buff,fasta]=sys.argv #assume the varfile has *.bed in the end
# Usage
# python makeigvpesr_trio.py varfile fasta sample ped cram_list buffer chromosome
# bash IL.DUP.HG00514.V2.sh
# bash igv.sh -b IL.DUP.HG00514.V2.txt

parser = argparse.ArgumentParser("makeigvsplit_trio.py")
parser.add_argument('-v', '--varfile', type=str, help='variant file including CHR, POS, END and SVID')
parser.add_argument('-fam_id','--fam_id', type=str, help='family to plot')
parser.add_argument('-p', '--ped', type=str, help='ped file')
#parser.add_argument('cram_list', type=str, help='comma separated list of all cram files to run igv on')
parser.add_argument('-samples', '--samples', type=str, help='List of all samples to run igv on')
parser.add_argument('-pe', '--pe', type=str, help='File of all pe files to run igv on')
parser.add_argument('-sr', '--sr', type=str, help='File of all sr files to run igv on')
parser.add_argument('-o', '--outdir', type=str, help = 'output folder')
parser.add_argument('-b', '--buff', type=str, help='length of buffer to add around variants', default=500)
parser.add_argument('-l', '--large_buff', type=str, help='length of buffer for large regions to add around variants', default=500)
parser.add_argument('-c', '--chromosome', type=str, help='name of chromosome to make igv on', default='all')
parser.add_argument('-i', '--igvfile', type=str, help='name of chromosome to make igv on', default='all')
parser.add_argument('-bam', '--bamfiscript', type=str, help='name of chromosome to make igv on', default='all')

args = parser.parse_args()

buff = int(args.buff)
large_buff = int(args.large_buff)
#fasta = args.fasta
varfile = args.varfile
pedigree = args.ped
fam_id = args.fam_id

outstring=os.path.basename(varfile)[0:-4]
bamdir="pe_bam"
outdir=args.outdir
igvfile=args.igvfile
bamfiscript=args.bamfiscript
###################################

#crams = args.crams
chromosome = args.chromosome

def ped_info_readin(ped_file):
    out={}
    fin=open(ped_file)
    for line in fin:
        pin=line.strip().split()
        if not pin[1] in out.keys():
            out[pin[1]]=[pin[1]]
        if not(pin[2])==0:
            out[pin[1]].append(pin[2])
        if not(pin[3])==0:
            out[pin[1]].append(pin[3])
    fin.close()
    return out

def cram_info_readin(cram_file):
    out={}
    fin=open(cram_file)
    for line in fin:
        pin=line.strip().split()
        if not pin[0] in out.keys():
            out[pin[0]]=pin[1:]
    fin.close()
    return(out)

#ped_info = ped_info_readin(args.ped)
#cram_info = cram_info_readin(args.cram_list)

#If file inputs
pe_colnames = colnames=[ 'pe'] 
pe = pd.read_csv(args.pe, sep='\t', names= pe_colnames, header=None).replace(np.nan, '', regex=True)
pe_list = pe['pe'].tolist()

sr_colnames = colnames=[ 'sr'] 
sr = pd.read_csv(args.sr, sep='\t', names= sr_colnames, header=None).replace(np.nan, '', regex=True)
sr_list = sr['sr'].tolist()
#cram_list = [c.replace('gs://', '/cromwell_root/') for c in cram_list_]

#sample_colnames = colnames=[ 'samples'] 
#sample = pd.read_csv(args.samples, sep='\t', names= sample_colnames, header=None).replace(np.nan, '', regex=True)
#samples_list = sample['samples'].tolist()

samples_list = args.samples.split(',')
#cram_list=args.crams.split(',')
'''
mydict = {key:value for key, value in zip(samples_list,cram_list)}
ped = pd.read_csv(pedigree, sep='\t', header=0).replace(np.nan, '', regex=True)
ped['FatherID'] = ped['FatherID'].astype(str)
ped['MotherID'] = ped['MotherID'].astype(str)
ped.Affected = pd.to_numeric(ped.Affected)
for sample_id in samples_list:
	if(ped.loc[(ped['IndividualID'] == sample_id)]['Affected'].iloc[0] == 2):
		if((ped.loc[(ped['IndividualID'] == sample_id)]['MotherID'].iloc[0] != '0' )| (ped.loc[(ped['IndividualID'] == sample_id)]['FatherID'].iloc[0] != '0' )):
			print(sample_id)
			proband_cram_file = mydict[sample_id]
			cram_list.remove(proband_cram_file)
			cram_list.insert(0, proband_cram_file)
		else:
			affected_cram_file = mydict[sample_id]
			cram_list.remove(affected_cram_file)
			cram_list.insert(1, affected_cram_file)
print(cram_list)
'''

with open(bamfiscript,'w') as h:
    h.write("#!/bin/bash\n")
    h.write("set -e\n")
    h.write("mkdir -p {}\n".format(bamdir))
    h.write("mkdir -p {}\n".format(outdir))
    with open(igvfile,'w') as g:
        g.write('new\n')
        with open(varfile,'r') as f:
            for line in f:
                dat=line.rstrip().split("\t")
                Chr=dat[0]
                if not chromosome=='all':
                    if not Chr == chromosome: continue
                Start_Buff=str(int(dat[1])-buff)
                End_Buff=str(int(dat[2])+buff)
                Start=str(int(dat[1]))
                End=str(int(dat[2]))
                ID=dat[3]
                for pe in pe_list:
                        g.write('load '+pe+'\n')
                for sr in sr_list:
                        g.write('load '+sr+'\n')
                if int(End)-int(Start)<10000:
                    g.write('goto '+Chr+":"+Start_Buff+'-'+End_Buff+'\n')
                    g.write('region '+Chr+":"+Start+'-'+End+'\n')
                    g.write('sort base\n')
                    g.write('viewaspairs\n')
                    #g.write('squish\n')
                    g.write('collapse Refseq Genes\n')
                    g.write('squish '+pe+'\n')
                    g.write('squish '+sr+'\n')
                    g.write('snapshotDirectory '+outdir+'\n')
                    g.write('snapshot '+fam_id+'_'+ID+'.png\n' )
                else:
                    g.write('goto '+Chr+":"+Start_Buff+'-'+str(int(Start_Buff)+large_buff)+'\n') # Extra 1kb buffer if variant large
                    g.write('region '+Chr+":"+Start+'-'+str(int(Start))+'\n') 
                    g.write('sort base\n')
                    g.write('viewaspairs\n')
                    #g.write('squish\n')
                    g.write('collapse Refseq Genes\n')
                    g.write('squish '+pe+'\n')
                    g.write('squish '+sr+'\n')
                    g.write('snapshotDirectory '+outdir+'\n')
                    g.write('snapshot '+fam_id+'_'+ID+'.left.png\n' )
                    g.write('goto '+Chr+":"+str(int(End)-large_buff)+'-'+End_Buff+'\n')
                    g.write('region '+Chr+":"+str(int(End))+'-'+End+'\n')
                    g.write('sort base\n')
                    g.write('viewaspairs\n')
                    #g.write('squish\n')
                    g.write('collapse Refseq Genes\n')
                    g.write('squish '+pe+'\n')
                    g.write('squish '+sr+'\n')
                    g.write('snapshotDirectory '+outdir+'\n')
                    g.write('snapshot '+fam_id+'_'+ID+'.right.png\n' )
                g.write('new\n')
        g.write('exit\n')