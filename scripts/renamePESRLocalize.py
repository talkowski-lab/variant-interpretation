import pandas as pd
import argparse
from subprocess import Popen, PIPE
import subprocess
import os



parser = argparse.ArgumentParser(description='Parse arguments')
parser.add_argument('--ped', dest='ped', help='Output file')
parser.add_argument('--pesr', dest='pesr', help='Output file')

args = parser.parse_args()

ped_file = args.ped
sample_cc = args.pesr

'''
def mv(file, output):
	#process1 = Popen(['gcloud', 'auth', 'application-default', 'print-access-token'], stdout=PIPE)
	#os.environ['GCS_OAUTH_TOKEN'] = process1.stdout.read().decode()
	result = subprocess.run(["mv", file, output])
	print(result.stderr)
'''

#find number of trios in ped file
ped=pd.read_csv(ped_file,sep='\t')
colnames = colnames=['sample', 'pe', 'sr']
pe_sr=pd.read_csv(sample_cc,sep='\t',names=colnames)
pd.options.display.max_colwidth = 500

#change to cromwell root so they can be moved
pe_sr['pe'] = pe_sr['pe'].str.replace('gs://', '/cromwell_root/')
pe_sr['sr'] = pe_sr['sr'].str.replace('gs://', '/cromwell_root/')

'''
cram = crams
for col in crams.columns:
	cram[col] = crams[col].apply(lambda x: x.strip())
print(cram['cram'])
'''
maternal_list = ped["MotherID"]
paternal_list = ped["FatherID"]
sample_list = ped["IndividualID"]

maternal_list = [i for i in maternal_list if i != '0']
print("Number of moms: ",len(maternal_list))
sr_file_mom = pe_sr.loc[pe_sr['sample'].isin(maternal_list), 'sr'].to_string(index = False)
pe_file_mom = pe_sr.loc[pe_sr['sample'].isin(maternal_list), 'pe'].to_string(index = False)
sr_file_mom_suffix = sr_file_mom.split('/')[-1]
pe_file_mom_suffix = pe_file_mom.split('/')[-1]
sr_file_mom_prefix = sr_file_mom.split(sr_file_mom_suffix)[0]
pe_file_mom_prefix = pe_file_mom.split(pe_file_mom_suffix)[0]
#mv(cram_file_mom, cram_file_mom_prefix + 'MOTHER.' + str(cram_file_mom_suffix))
#mv(crai_file_mom, crai_file_mom_prefix +  'MOTHER.' + str(crai_file_mom_suffix))
pe_sr.loc[pe_sr['sample'].isin(maternal_list), 'new_sr'] = 'MOTHER.' + str(sr_file_mom_suffix)
pe_sr.loc[pe_sr['sample'].isin(maternal_list), 'new_pe'] = 'MOTHER.' + str(pe_file_mom_suffix)


paternal_list = [i for i in paternal_list if i != '0']
print("Number of dads: ", len(paternal_list))
sr_file_dad = pe_sr.loc[pe_sr['sample'].isin(paternal_list), 'sr'].to_string(index = False)
pe_file_dad = pe_sr.loc[pe_sr['sample'].isin(paternal_list), 'pe'].to_string(index = False)
sr_file_dad_suffix = sr_file_dad.split('/')[-1]
pe_file_dad_suffix = pe_file_dad.split('/')[-1]
sr_file_dad_prefix = sr_file_dad.split(sr_file_dad_suffix)[0]
pe_file_dad_prefix = pe_file_dad.split(pe_file_dad_suffix)[0]
#mv(cram_file_dad, cram_file_dad_prefix + 'FATHER.' + str(cram_file_dad_suffix))
#mv(crai_file_dad, crai_file_dad_prefix + 'FATHER.' + str(crai_file_dad_suffix))
pe_sr.loc[pe_sr['sample'].isin(paternal_list), 'new_sr'] = 'FATHER.' + str(sr_file_dad_suffix)
pe_sr.loc[pe_sr['sample'].isin(paternal_list), 'new_pe'] = 'FATHER.' + str(pe_file_dad_suffix)


proband_list = [i for i in sample_list if (i not in maternal_list) & (i not in paternal_list)]
print("Number of probands: ", len(proband_list))



#proband = ped[ped['subject_id'].isin([proband_list])]
mask = ped.IndividualID.astype(str).isin(proband_list)
proband = ped[mask]


proband.Affected = pd.to_numeric(proband.Affected)
affected_proband = proband[proband['Affected'] == 2]

affected_list_proband = affected_proband['IndividualID']
print("Number of affected probands: ",len(affected_list_proband))
sr_file_p = pe_sr.loc[pe_sr['sample'].isin(affected_list_proband), 'sr'].to_string(index = False)
pe_file_p = pe_sr.loc[pe_sr['sample'].isin(affected_list_proband), 'pe'].to_string(index = False)
sr_file_p_suffix = sr_file_p.split('/')[-1]
pe_file_p_suffix = pe_file_p.split('/')[-1]
sr_file_p_prefix = sr_file_p.split(sr_file_p_suffix)[0]
pe_file_p_prefix = pe_file_p.split(pe_file_p_suffix)[0]
#mv(cram_file_p, cram_file_p_prefix + 'PROBAND.' + str(cram_file_p_suffix))
#mv(crai_file_p, crai_file_p_prefix + 'PROBAND.' + str(crai_file_p_suffix))
pe_sr.loc[pe_sr['sample'].isin(affected_list_proband), 'new_sr'] = 'PROBAND.' + str(sr_file_p_suffix)
pe_sr.loc[pe_sr['sample'].isin(affected_list_proband), 'new_pe'] = 'PROBAND.' + str(pe_file_p_suffix)


sibling = proband[proband['Affected'] == 1]
sibling_list = sibling['IndividualID']
print("Number of siblings: ",len(sibling_list))
sr_file_sib = pe_sr.loc[pe_sr['sample'].isin(sibling_list), 'sr'].to_string(index = False)
pe_file_sib = pe_sr.loc[pe_sr['sample'].isin(sibling_list), 'pe'].to_string(index = False)
sr_file_sib_suffix = sr_file_sib.split('/')[-1]
pe_file_sib_suffix = pe_file_sib.split('/')[-1]
sr_file_sib_prefix = sr_file_sib.split(sr_file_sib_suffix)[0]
pe_file_sib_prefix = pe_file_sib.split(pe_file_sib_suffix)[0]
#mv(cram_file_sib, cram_file_sib_prefix + 'SIBLING.' + str(cram_file_sib_suffix))
#mv(crai_file_sib, cram_file_sib_prefix + 'SIBLING.' + str(crai_file_sib_suffix))
pe_sr.loc[pe_sr['sample'].isin(sibling_list), 'new_sr'] = 'SIBLING.' + str(sr_file_sib_suffix)
pe_sr.loc[pe_sr['sample'].isin(sibling_list), 'new_pe'] = 'SIBLING.' + str(pe_file_sib_suffix)


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
pe_sr.to_csv(path_or_buf='changed_sample_pe_sr.txt', mode='a', index=False, sep='\t', header=False)