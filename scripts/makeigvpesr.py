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
#parser.add_argument('-n', '--nestedrepeats', type=str, help='nested repeats sequences')
#parser.add_argument('-s', '--simplerepeats', type=str, help='simple repeats sequences')
#parser.add_argument('-e', '--emptytrack', type=str, help='empty track')
#parser.add_argument('-f', '--fasta', type=str, help='reference sequences')
#parser.add_argument('sample', type=str, help='name of sample to make igv on')
parser.add_argument('-fam_id','--fam_id', type=str, help='family to plot')
parser.add_argument('-p', '--ped', type=str, help='ped file')
#parser.add_argument('cram_list', type=str, help='comma separated list of all cram files to run igv on')
parser.add_argument('-samples', '--samples', type=str, help='List of all samples to run igv on')
parser.add_argument('-crams', '--crams', type=str, help='File of all cram files to run igv on')
parser.add_argument('-o', '--outdir', type=str, help = 'output folder')
parser.add_argument('-b', '--buff', type=str, help='length of buffer to add around variants', default=500)
parser.add_argument('-c', '--chromosome', type=str, help='name of chromosome to make igv on', default='all')
parser.add_argument('-i', '--igvfile', type=str, help='name of chromosome to make igv on', default='all')
parser.add_argument('-bam', '--bamfiscript', type=str, help='name of chromosome to make igv on', default='all')
parser.add_argument('-m', '--igvmaxwindow', type=str, help='max length of SV to appear in IGV', default=10e10)

args = parser.parse_args()


buff = int(args.buff)
#fasta = args.fasta
varfile = args.varfile
pedigree = args.ped
fam_id = args.fam_id
igv_max_window = args.igvmaxwindow


outstring=os.path.basename(varfile)[0:-4]
bamdir="pe_bam"
outdir=args.outdir
igvfile=args.igvfile
bamfiscript=args.bamfiscript
###################################

#crams = args.crams
chromosome = args.chromosome
#nested_repeats = args.nestedrepeats
#simple_repeats = args.simplerepeats
#empty_track = args.emptytrack

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
cram_colnames = colnames=[ 'cram']
cram = pd.read_csv(args.crams, sep='\t', names= cram_colnames, header=None).replace(np.nan, '', regex=True)
cram_list = cram['cram'].tolist()
#cram_list = [c.replace('gs://', '/cromwell_root/') for c in cram_list_]

#sample_colnames = colnames=[ 'samples']
#sample = pd.read_csv(args.samples, sep='\t', names= sample_colnames, header=None).replace(np.nan, '', regex=True)
#samples_list = sample['samples'].tolist()

samples_list = args.samples.split(',')
#cram_list=args.crams.split(',')
mydict = {key:value for key, value in zip(samples_list,cram_list)}
ped=pd.read_csv(pedigree,sep='\t').iloc[:,:6]
ped.columns = ['FamilyID','IndividualID','MotherID','FatherID','Sex','Affected']
ped[['FamilyID','IndividualID','MotherID','FatherID']] = ped[['FamilyID','IndividualID','MotherID','FatherID']].astype(str)

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
                Start=int(dat[1])
                End=int(dat[2])
                ID=dat[3]
                Length=End-Start

                Length_total=int(Length+(Length)*1.5)

                for cram in cram_list:
                        g.write('load '+cram+'\n')

                if Length_total<int(igv_max_window):
                    if Length_total<1000:
                        Start_Buff=int(Start-500)
                        End_Buff=int(End+500)
                    else:
                        Start_Buff = int(Start - (Length * 0.25))
                        End_Buff = int(End + (Length * 0.25))
                    g.write('goto '+Chr+":"+str(Start_Buff)+'-'+str(End_Buff)+'\n')
                    g.write('region '+Chr+":"+str(Start)+'-'+str(End)+'\n')
                    g.write('sort base\n')
                    g.write('viewaspairs\n')
                    g.write('squish\n')
                    g.write('collapse Refseq Genes\n')
                    g.write('snapshotDirectory '+outdir+'\n')
                    g.write('snapshot '+fam_id+'_'+ID+'.png\n' )
                else:
                    g.write('goto '+Chr+":"+str(Start-buff)+'-'+str(Start+buff)+'\n')
                    g.write('region '+Chr+":"+str(Start)+'-'+str(Start)+'\n')
                    g.write('sort base\n')
                    g.write('viewaspairs\n')
                    g.write('squish\n')
                    g.write('collapse Refseq Genes\n')
                    g.write('snapshotDirectory '+outdir+'\n')
                    g.write('snapshot '+fam_id+'_'+ID+'.left.png\n' )
                    g.write('goto '+Chr+":"+str(End-buff)+'-'+str(End+buff)+'\n')
                    g.write('region '+Chr+":"+str(End)+'-'+str(End)+'\n')
                    g.write('sort base\n')
                    g.write('viewaspairs\n')
                    g.write('squish\n')
                    g.write('collapse Refseq Genes\n')
                    g.write('snapshotDirectory '+outdir+'\n')
                    g.write('snapshot '+fam_id+'_'+ID+'.right.png\n' )
                # g.write('goto '+Chr+":"+Start+'-'+End+'\n')
                # g.write('sort base\n')
                # g.write('viewaspairs\n')
                # g.write('squish\n')
                # g.write('snapshotDirectory '+outdir+'\n')
                # g.write('snapshot '+ID+'.png\n' )
                g.write('new\n')
        g.write('exit\n')