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


#def mv(file, output):
	#process1 = Popen(['gcloud', 'auth', 'application-default', 'print-access-token'], stdout=PIPE)
	#os.environ['GCS_OAUTH_TOKEN'] = process1.stdout.read().decode()
	#result = subprocess.run(["mv", file, output])
	#print(result.stderr)

#find number of trios in ped file
ped=pd.read_csv(ped_file,sep='\t')
ped[['FamilyID','IndividualID','MotherID','FatherID']] = ped[['FamilyID','IndividualID','MotherID','FatherID']].astype(str)

colnames = colnames=['sample', 'crai', 'cram']
cram=pd.read_csv(sample_cc,sep='\t',names=colnames)
pd.options.display.max_colwidth = 500

#change to cromwell root so they can be moved
#cram['cram'] = cram['cram'].str.replace('gs://', '/cromwell_root/')
#cram['crai'] = cram['crai'].str.replace('gs://', '/cromwell_root/')


maternal_list = ped["MotherID"]
paternal_list = ped["FatherID"]
sample_list = ped["IndividualID"]

maternal_list = [i for i in maternal_list if i != '0']
print("Number of moms: ",len(maternal_list))
cram_file_mom = cram.loc[cram['sample'].isin(maternal_list), 'cram'].to_string(index = False)
crai_file_mom = cram.loc[cram['sample'].isin(maternal_list), 'crai'].to_string(index = False)
cram_file_mom_suffix = cram_file_mom.split('/')[-1]
crai_file_mom_suffix = crai_file_mom.split('/')[-1]
cram_file_mom_prefix = cram_file_mom.split(cram_file_mom_suffix)[0]
crai_file_mom_prefix = crai_file_mom.split(crai_file_mom_suffix)[0]
#mv(cram_file_mom, 'mother.cram')
#mv(crai_file_mom, 'mother.cram.crai')
cram.loc[cram['sample'].isin(maternal_list), 'new_cram'] = cram_file_mom_prefix + 'MOTHER.' + str(cram_file_mom_suffix)
cram.loc[cram['sample'].isin(maternal_list), 'new_crai'] = crai_file_mom_prefix + 'MOTHER.' + str(crai_file_mom_suffix)


paternal_list = [i for i in paternal_list if i != '0']
print("Number of dads: ", len(paternal_list))
cram_file_dad = cram.loc[cram['sample'].isin(paternal_list), 'cram'].to_string(index = False)
crai_file_dad = cram.loc[cram['sample'].isin(paternal_list), 'crai'].to_string(index = False)
cram_file_dad_suffix = cram_file_dad.split('/')[-1]
crai_file_dad_suffix = crai_file_dad.split('/')[-1]
cram_file_dad_prefix = cram_file_dad.split(cram_file_dad_suffix)[0]
crai_file_dad_prefix = crai_file_dad.split(crai_file_dad_suffix)[0]
#mv(cram_file_dad, 'father.cram')
#mv(crai_file_dad, 'father.cram.crai')
cram.loc[cram['sample'].isin(paternal_list), 'new_cram'] = cram_file_dad_prefix + 'FATHER.' + str(cram_file_dad_suffix)
cram.loc[cram['sample'].isin(paternal_list), 'new_crai'] = crai_file_dad_prefix + 'FATHER.' + str(crai_file_dad_suffix)


proband_list = [i for i in sample_list if (i not in maternal_list) & (i not in paternal_list)]
print("Number of probands: ", len(proband_list))



#proband = ped[ped['subject_id'].isin([proband_list])]
mask = ped.IndividualID.astype(str).isin(proband_list)
proband = ped[mask]


proband.Affected = pd.to_numeric(proband.Affected)
affected_proband = proband[proband['Affected'] == 2]

affected_list_proband = affected_proband['IndividualID']
print("Number of affected probands: ",len(affected_list_proband))
cram_file_p = cram.loc[cram['sample'].isin(affected_list_proband), 'cram'].to_string(index = False)
crai_file_p = cram.loc[cram['sample'].isin(affected_list_proband), 'crai'].to_string(index = False)
cram_file_p_suffix = cram_file_p.split('/')[-1]
crai_file_p_suffix = crai_file_p.split('/')[-1]
cram_file_p_prefix = cram_file_p.split(cram_file_p_suffix)[0]
crai_file_p_prefix = crai_file_p.split(crai_file_p_suffix)[0]
#mv(cram_file_p, 'proband.' + str(cram_file_p))
#mv(crai_file_p, 'proband.' + str(crai_file_p))
cram.loc[cram['sample'].isin(affected_list_proband), 'new_cram'] = cram_file_p_prefix + 'PROBAND.' + str(cram_file_p_suffix)
cram.loc[cram['sample'].isin(affected_list_proband), 'new_crai'] = crai_file_p_prefix + 'PROBAND.' + str(crai_file_p_suffix)


sibling = proband[proband['Affected'] == 1]
sibling_list = sibling['IndividualID']
print("Number of siblings: ",len(sibling_list))
cram_file_sib = cram.loc[cram['sample'].isin(sibling_list), 'cram'].to_string(index = False)
crai_file_sib = cram.loc[cram['sample'].isin(sibling_list), 'crai'].to_string(index = False)
cram_file_sib_suffix = cram_file_sib.split('/')[-1]
crai_file_sib_suffix = crai_file_sib.split('/')[-1]
cram_file_sib_prefix = cram_file_sib.split(cram_file_sib_suffix)[0]
crai_file_sib_prefix = crai_file_sib.split(crai_file_sib_suffix)[0]
#mv(cram_file_sib, 'sibling.cram')
#mv(crai_file_sib, 'sibling.cram.crai')
cram.loc[cram['sample'].isin(sibling_list), 'new_cram'] = cram_file_sib_prefix + 'SIBLING.' + str(cram_file_sib_suffix)
cram.loc[cram['sample'].isin(sibling_list), 'new_crai'] = cram_file_sib_prefix + 'SIBLING.' + str(crai_file_sib_suffix)


mat = {

	'MotherID': maternal_list,
}

mat_df = pd.DataFrame(mat)

pat = {

	'FatherID': paternal_list,
}

pat_df = pd.DataFrame(pat)

sib = {

	'SiblingID': sibling_list,
}

sib_df = pd.DataFrame(sib)

affected_list = {

	'SiblingID': affected_list_proband,
}

affected_df = pd.DataFrame(affected_list)


#mat_df.to_csv(path_or_buf='mothers.txt', mode='a', index=False, sep='\t', header=False)
#pat_df.to_csv(path_or_buf='fathers.txt', mode='a', index=False, sep='\t', header=False)
#sib_df.to_csv(path_or_buf='siblings.txt', mode='a', index=False, sep='\t', header=False)
#affected_df.to_csv(path_or_buf='affected_proband.txt', mode='a', index=False, sep='\t', header=False)
cram.to_csv(path_or_buf='changed_sample_crai_cram.txt', mode='a', index=False, sep='\t', header=False)