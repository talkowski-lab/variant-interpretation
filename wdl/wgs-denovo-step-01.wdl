version 1.0

workflow step1 {
	input {
		File python_trio_sample_script
		File python_preprocess_script
		File lcr_uri
		File ped_uri
		Array[Array[File]] vcf_uri_list
		File? meta_uri
		File? trio_uri
		String cohort_prefix
		String hail_docker
	}

	if (!defined(meta_uri)) {
		call makeTrioSampleFilesNew {
			input:
				python_trio_sample_script=python_trio_sample_script,
				ped_uri=ped_uri,
				cohort_prefix=cohort_prefix,
				hail_docker=hail_docker
		}
		File meta_uri = select_first([makeTrioSampleFiles.meta_uri_out])
		File trio_uri = select_first([makeTrioSampleFiles.trio_uri_out])
	}

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

task makeTrioSampleFilesNew {
	input {
		File python_trio_sample_script
		File ped_uri
		String cohort_prefix
		String hail_docker
	}

	runtime {
		docker: hail_docker
	}

	output {
		File meta_uri_out = "${cohort_prefix}_sample_list.txt"
		File trio_uri_out = "${cohort_prefix}_trio_list.txt"
	}

	command <<<
	python3 ${python_trio_sample_script} ${ped_uri} ${cohort_prefix}
	>>>
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

	output {
		File preprocessed_vcf = basename(vcf_uri, '.vcf.gz') + '.filtered.vcf.gz'
	}

	command <<<
	python3 ${python_preprocess_script} ${lcr_uri} ${ped_uri} ${meta_uri} ${trio_uri} ${vcf_uri}
	>>>
}
