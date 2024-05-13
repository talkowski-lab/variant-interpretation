import pandas as pd
import argparse
from subprocess import Popen, PIPE
import subprocess
import os



parser = argparse.ArgumentParser(description='Parse arguments')
parser.add_argument('--ped', dest='ped', help='Output file')
parser.add_argument('--scc', dest='scc', help='Output file')

args = parser.parse_args()

ped_file = args.ped
sample_cc = args.scc

'''
def mv(file, output):
	#process1 = Popen(['gcloud', 'auth', 'application-default', 'print-access-token'], stdout=PIPE)
	#os.environ['GCS_OAUTH_TOKEN'] = process1.stdout.read().decode()
	result = subprocess.run(["mv", file, output])
	print(result.stderr)
'''

#find number of trios in ped file
#find number of trios in ped file
ped=pd.read_csv(ped_file,sep='\t').iloc[:,:6]
ped.columns = ['FamilyID','IndividualID','MotherID','FatherID','Sex','Affected']
ped[['FamilyID','IndividualID','MotherID','FatherID']] = ped[['FamilyID','IndividualID','MotherID','FatherID']].astype(str)

colnames = colnames=['sample', 'crai', 'cram']
cram=pd.read_csv(sample_cc,sep='\t',names=colnames)
cram[colnames] = cram[colnames].astype(str)

fathers = ped[~ped.FatherID.isin(['0','-9'])].set_index('IndividualID').FatherID.to_dict()
mothers = ped[~ped.MotherID.isin(['0','-9'])].set_index('IndividualID').MotherID.to_dict()

def get_sample_role(row):
    if (row.MotherID=='0') & (row.FatherID=='0'):
        if (row.IndividualID in fathers.values()):
            role = 'FATHER'
        elif (row.IndividualID in mothers.values()):
            role = 'MOTHER'
        else:
            role = 'UNKNOWN'
    elif (row.MotherID in mothers.values()) & (row.FatherID in fathers.values()):
        role = 'PROBAND'
    else:
        role = 'SIBLING'
    return role

ped['role'] = ped.apply(get_sample_role, axis=1)
cram['role'] = cram['sample'].map(ped.set_index('IndividualID').role.to_dict())

#change to cromwell root so they can be moved
cram['cram'] = cram['cram'].str.replace('gs://', '/cromwell_root/')
cram['crai'] = cram['crai'].str.replace('gs://', '/cromwell_root/')

cram['new_cram'] = cram.cram.apply(os.path.dirname) + '/' + cram.role + '.' + cram.cram.apply(os.path.basename)
cram['new_crai'] = cram.crai.apply(os.path.dirname) + '/' + cram.role + '.' + cram.crai.apply(os.path.basename)

cram = cram.drop('role', axis=1)
cram.to_csv(path_or_buf='changed_sample_crai_cram.txt', index=False, sep='\t', header=False)