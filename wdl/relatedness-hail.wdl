version 1.0

import "mergeVCFs.wdl" as mergeVCFs
import "wes-denovo-helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow Relatedness {
    input {
        Array[File] vep_vcf_files
        File? merged_vep_file
        File ped_uri
        File bed_file
        String cohort_prefix
        String sv_base_mini_docker
        String hail_docker
    }

    if (!defined(merged_vep_file)) {
        
        scatter (vcf_uri in vep_vcf_files) {
            String filename = basename(vcf_uri)
            String prefix = if (sub(filename, ".gz", "")!=filename) then basename(filename, ".vcf.gz") else basename(filename, ".vcf.bgz")
            call helpers.subsetVCFs as subsetVCFs {
                input:
                    bed_file=bed_file,
                    vcf_uri=vcf_uri,
                    vcf_idx=vcf_uri+'.tbi',
                    output_name=prefix + '.somalier.subset.vcf.gz',
                    sv_base_mini_docker=sv_base_mini_docker
            }
        }

        call mergeVCFs.mergeVCFs as mergeVCFs {
        input:
            vcf_files=subsetVCFs.subset_vcf,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix
        }
    }

    # File merged_vcf_file = select_first([merged_vep_file, mergeVCFs.merged_vcf_file])

    output {
    File merged_vcf_file = select_first([merged_vep_file, mergeVCFs.merged_vcf_file])
    }
}