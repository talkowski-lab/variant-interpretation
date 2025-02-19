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

workflow FilterRegion {
    input {
        Array[File] vcf_files
        File bed_file
        String sv_base_mini_docker
        String sort_after_merge
        String prefix
        RuntimeAttr? runtime_attr_filter_vcf
        RuntimeAttr? runtime_attr_merge_vcfs
    }

    scatter (vcf_uri in select_first([vcf_files])) {
        String filename = basename(vcf_uri)
        call helpers.subsetVCFs as subsetVCFs {
            input:
                bed_file=bed_file,
                vcf_uri=vcf_uri,
                vcf_idx=vcf_uri+'.tbi',
                output_name=prefix + '.subset.vcf.gz',
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_filter_vcf
        }
    }

    call mergeVCFs.mergeVCFs as mergeVCFs {
        input:
            vcf_files=subsetVCFs.subset_vcf,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=prefix,
            sort_after_merge=sort_after_merge,
            runtime_attr_override=runtime_attr_merge_vcfs
        }

    output {
        File merged_vcf_file = mergeVCFs.merged_vcf_file
    }
}