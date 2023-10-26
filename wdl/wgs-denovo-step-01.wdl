version 1.0

workflow step1 {
	input {
		File python_trio_sample_script
		File python_preprocess_script
		File lcr_uri
		File ped_uri
		Array[Array[File]] vcf_uri_list
		String bucket_id
		String cohort_prefix
		String hail_docker
	}

	call makeTrioSampleFiles {
		input:
			python_trio_sample_script=python_trio_sample_script,
			ped_uri=ped_uri,
			bucket_id=bucket_id,
			cohort_prefix=cohort_prefix,
			hail_docker=hail_docker
	}
	File meta_uri = makeTrioSampleFiles.meta_uri
	File trio_uri = makeTrioSampleFiles.trio_uri

	scatter (vcf_uri_sublist in vcf_uri_list) {
		scatter (vcf_uri in vcf_uri_sublist) {
			call preprocessVCF {
				input:
					python_preprocess_script=python_preprocess_script,
					lcr_uri=lcr_uri,
					ped_uri=ped_uri,
					vcf_uri=vcf_uri,
					meta_uri=meta_uri,
					trio_uri=trio_uri,
					hail_docker=hail_docker
			}
		}
	}

	output {
		Array[Array[File]] preprocessed_vcf_list = preprocessVCF.preprocessed_vcf
	}
}

task makeTrioSampleFiles {
	input {
		File python_trio_sample_script
		File ped_uri
		String bucket_id
		String cohort_prefix
		String hail_docker
	}

	runtime {
		docker: hail_docker
	}

	command <<<
	python3 ~{python_trio_sample_script} ~{ped_uri} ~{cohort_prefix}
	>>>
	
	output {
		File meta_uri = "~{bucket_id}/resources/metadata/~{cohort_prefix}_sample_list.txt"
		File trio_uri = "~{bucket_id}/resources/metadata/~{cohort_prefix}_trio_list.txt"
	}
}

task preprocessVCF {
	input {
		File python_preprocess_script
		File lcr_uri
		File ped_uri
		File vcf_uri
		File meta_uri
		File trio_uri
		String hail_docker
	}

	runtime {
		docker: hail_docker
	}

	command <<<
	python3 ~{python_preprocess_script} ~{lcr_uri} ~{ped_uri} ~{meta_uri} ~{trio_uri} ~{vcf_uri}
	>>>

	output {
		File preprocessed_vcf = basename(vcf_uri, '.vcf.gz') + '.filtered.vcf.gz'
	}
}
