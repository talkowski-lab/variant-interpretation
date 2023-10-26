version 1.0

workflow makeTrioSampleFiles {
	input {
		# File python_script
		File ped_uri
		String cohort_prefix
		String hail_docker
	}

	call makeTrioSampleFilesTask {
		input:
			ped_uri=ped_uri,
			cohort_prefix=cohort_prefix,
			hail_docker=hail_docker
	}

	output {
		File meta_file = "${cohort_prefix}_sample_list.txt"
		File trio_file = "${cohort_prefix}_trio_list.txt"
	}
}

task makeTrioSampleFilesTask {
	input {
		File python_script
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
	python3 ${python_script} ${ped_uri} ${cohort_prefix}
	>>>
}
