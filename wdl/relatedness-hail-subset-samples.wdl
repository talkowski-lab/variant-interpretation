version 1.0

import "mergeVCFs.wdl" as mergeVCFs
import "wes-denovo-helpers.wdl" as helpers
import "relatedness-hail.wdl" as relatednessHail

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
        Int samples_per_chunk
        String cohort_prefix
        String relatedness_qc_script
        String plot_relatedness_script
        String sex_qc_script
        String sv_base_mini_docker
        String hail_docker
        String bucket_id
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
    File merged_vcf_file = select_first([merged_vep_file, mergeVCFs.merged_vcf_file])

    call relatednessHail.imputeSex as imputeSex {
        input:
        vcf_uri=merged_vcf_file,
        bed_file=bed_file,
        ped_uri=ped_uri,
        cohort_prefix=cohort_prefix,
        sex_qc_script=sex_qc_script,
        hail_docker=hail_docker
    }

    call helpers.splitSamples as splitSamples {
        input:
            vcf_file=merged_vcf_file,
            samples_per_chunk=samples_per_chunk,
            cohort_prefix=cohort_prefix,
            sv_base_mini_docker=sv_base_mini_docker
    } 

    scatter (sample_file in splitSamples.sample_shard_files) {
        call helpers.subsetVCFSamplesHail as subsetVCFSamples {
            input:
                samples_file=sample_file,
                vcf_file=merged_vcf_file,
                hail_docker=hail_docker
        }

        call relatednessHail.checkRelatedness as checkRelatedness {
            input:
            vcf_uri=subsetVCFSamples.vcf_subset,
            bed_file=bed_file,
            ped_uri=ped_uri,
            cohort_prefix=cohort_prefix,
            relatedness_qc_script=relatedness_qc_script,
            hail_docker=hail_docker,
            bucket_id=bucket_id
        }
    }
    call helpers.getHailMTSizes as getRelatednessQCSizes {
            input:
            mt_uris=checkRelatedness.relatedness_qc,
            hail_docker=hail_docker
    }
    call helpers.getHailMTSizes as getAnnotRelPedSizes {
            input:
            mt_uris=checkRelatedness.annot_rel_ped,
            hail_docker=hail_docker
    }

    call helpers.mergeResultsPython as mergeRelatednessQC {
            input:
            input_size=getAnnotRelPedSizes.mt_size,
            tsvs=checkRelatedness.relatedness_qc,
            merged_filename=cohort_prefix + '_relatedness_qc',
            hail_docker=hail_docker
    }

    call helpers.mergeResultsPython as mergeAnnotRelPed {
            input:
            input_size=getAnnotRelPedSizes.mt_size,
            tsvs=checkRelatedness.annot_rel_ped,
            merged_filename=cohort_prefix + '_annot_relationships',
            hail_docker=hail_docker
    }

    call relatednessHail.plotRelatedness as plotRelatedness {
        input:
        annot_rel_ped=mergeAnnotRelPed.merged_tsv,
        ped_uri=ped_uri,
        cohort_prefix=cohort_prefix,
        plot_relatedness_script=plot_relatedness_script,
        hail_docker=hail_docker
    }

    output {
        File sex_qc_plots = imputeSex.sex_qc_plots
        File ped_sex_qc = imputeSex.ped_sex_qc
        File relatedness_qc = mergeRelatednessQC.merged_tsv
        File annot_rel_ped = mergeAnnotRelPed.merged_tsv
        File relatedness_plot = plotRelatedness.relatedness_plot
    }
}