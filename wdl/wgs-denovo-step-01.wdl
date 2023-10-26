version 1.0

workflow step1 {
	input {
		# File python_script
		File ped_uri
		File? meta_file
		File? trio_file
		String cohort_prefix
		String hail_docker
	}

	if (!defined(meta_file)) {
		call makeTrioSampleFiles {
			input:
				ped_uri=ped_uri,
				cohort_prefix=cohort_prefix,
				hail_docker=hail_docker
		}
		File meta_file = makeTrioSampleFiles.meta_file
		File trio_file = makeTrioSampleFiles.trio_file
	}
	File meta_file = select_first[meta_file]
	File trio_file = select_first[trio_file]

	output {
	}
}

task makeTrioSampleFiles {
	input {
		File ped_uri
		String cohort_prefix
		String hail_docker
	}

	runtime {
		docker: hail_docker
	}

	output {
		File meta_file = "${cohort_prefix}_sample_list.txt"
		File trio_file = "${cohort_prefix}_trio_list.txt"
	}

	command <<<
	# python3 ${python_script} ${ped_uri}
	python3 <<CODE
	import pandas as pd
	import gcsfs
	import os
	import sys
	from google.cloud import storage

	ped_uri = "~{ped_uri}"
	cohort_prefix = "~{cohort_prefix}"
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
}
