version 1.0

workflow makeTrioSampleFiles {
	input {
		File python_script
		File ped_uri
		String cohort_prefix
		String hail_docker
	}

	call makeTrioSampleFilesTask {
		input:
			python_script=python_script,
			ped_uri=ped_uri,
			cohort_prefix=cohort_prefix,
			hail_docker=hail_docker
	}

	output {
		File meta_file = makeTrioSampleFilesTask.meta_file
		File trio_file = makeTrioSampleFilesTask.trio_file
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
		File meta_file = meta_file
		File trio_file = trio_file
	}

	command <<<
	python3 ~{python_script} ~{ped_uri} ~{cohort_prefix} --out1=trio_file --out2=meta_file
	>>>
}
