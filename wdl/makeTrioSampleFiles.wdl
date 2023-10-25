version 1.0

workflow makeTrioSampleFiles {
	input {
		# File python_script
		File ped_uri
		String cohort_prefix
		String hail_docker
	}

	output {
		meta_file = File "${cohort_prefix}_sample_list.txt"
		trio_file = File "${cohort_prefix}_trio_list.txt"
	}

	runtime {
		docker: hail_docker
	}

	command <<<
	# python3 ${python_script} ${ped_uri}
	python3 <<CODE
		import pandas as pd
		import gcsfs
		import os
		import sys
		from google.cloud import storage

		fs = gcsfs.GCSFileSystem(project='talkowski-sv-gnomad')
		with fs.open(ped_uri) as f:
	    	ped = pd.read_csv(f, sep='\t')
	    trio = ped[(ped.FatherID != '0') & (ped.MotherID != '0')].iloc[:, :4]
		trio['TrioID'] = trio['FamID'] + '-' + trio['IndividualID']
		trio.rename(columns={'IndividualID': 'SampleID'}, inplace=True)
		trio.to_csv(f"{cohort_prefix}_trio_list.txt", sep='\t', index=False)

		trio.columns.name = 'Role'
		trio.index = trio.FamID

		sample_data = trio.iloc[:, 1:-1].stack().reset_index()
		sample_data.columns = ['FamID', 'Role', 'SampleID']
		sample_data = sample_data.replace({'SampleID': 'child', 'MotherID': 'mother', 'FatherID': 'father'})
		sample_data.to_csv(f"{cohort_prefix}_sample_list.txt", sep='\t', index=False)

		CODE
		>>>
	>>>
}
