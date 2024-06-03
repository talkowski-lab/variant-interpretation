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

workflow runSomalier {
    input {
        File sites_uri
        File ref_fasta
        Array[File]? vep_vcf_files
        Array[File]? vep_annotated_final_vcf
        Array[String]? mt_shards
        File? merged_vcf_file
        File ped_uri
        File bed_file
        File ancestry_labels_1kg
        File somalier_1kg_tar
        String correct_somalier_ped_python_script
        String cohort_prefix
        String somalier_docker
        String sv_base_mini_docker
        String hail_docker
        Int samples_per_chunk
        Boolean subset_ped=true
        Boolean infer_ped=true
        Boolean unknown_flag=false
        Float relatedness_cutoff=0.25
        RuntimeAttr? runtime_attr_relatedness
        RuntimeAttr? runtime_attr_correct
    }

    if (!defined(merged_vcf_file)) {
        if (defined(mt_shards)) {
            scatter (mt_uri in select_first([mt_shards])) {
                call helpers.getHailMTSize {
                    input:
                        mt_uri=mt_uri,
                        hail_docker=hail_docker
                }
                call helpers.filterIntervalsToVCF as filterIntervals {
                    input:
                        bed_file=bed_file,
                        input_size=getHailMTSize.mt_size,
                        mt_uri=mt_uri,
                        hail_docker=hail_docker
                }
            }
        }

        if (!defined(mt_shards)) {
            Array[File] vep_files = select_first([vep_vcf_files, vep_annotated_final_vcf])
            
            scatter (vcf_uri in vep_files) {
                call subsetVCFs {
                    input:
                        bed_file=bed_file,
                        vcf_uri=vcf_uri,
                        vcf_idx=vcf_uri+'.tbi',
                        somalier_docker=somalier_docker
                }
            }
        }

        call mergeVCFs.mergeVCFs as mergeVCFs {
            input:
                vcf_files=select_first([filterIntervals.vcf_filt, subsetVCFs.subset_vcf]),
                sv_base_mini_docker=sv_base_mini_docker,
                cohort_prefix=cohort_prefix
        }
    }

    File merged_vcf_file_ = select_first([merged_vcf_file, mergeVCFs.merged_vcf_file])

    call helpers.splitSamples as splitSamples {
        input:
            vcf_file=merged_vcf_file_,
            samples_per_chunk=samples_per_chunk,
            cohort_prefix=cohort_prefix,
            sv_base_mini_docker=sv_base_mini_docker
    } 

    scatter (sample_file in splitSamples.sample_shard_files) {
        call helpers.subsetVCFSamplesHail {
            input:
                samples_file=sample_file,
                vcf_file=merged_vcf_file_,
                hail_docker=hail_docker
        }

        call relatedness_subset {
            input:
                sites_uri=sites_uri,
                ref_fasta=ref_fasta,
                vcf_uri=subsetVCFSamplesHail.vcf_subset,
                ped_uri=ped_uri,
                sample_file=sample_file,
                ancestry_labels_1kg=ancestry_labels_1kg,
                somalier_1kg_tar=somalier_1kg_tar,
                cohort_prefix=cohort_prefix,
                somalier_docker=somalier_docker,
                infer_ped=infer_ped,
                unknown_flag=unknown_flag,
                relatedness_cutoff=relatedness_cutoff,
                runtime_attr_override=runtime_attr_relatedness
        }
    }

    output {
        Array[File] out_samples = relatedness_subset.out_samples
        Array[File] out_pairs = relatedness_subset.out_pairs
        Array[File] out_groups = relatedness_subset.out_groups
        Array[File] out_html = relatedness_subset.out_html
        Array[File] ancestry_html = relatedness_subset.ancestry_html
        Array[File] ancestry_out = relatedness_subset.ancestry_out
    }
}

task subsetVCFs {
    input {
        File vcf_uri
        File vcf_idx
        File bed_file
        String somalier_docker
        RuntimeAttr? runtime_attr_override
    }

    Float relatedness_size = size(vcf_uri, "GB") 
    Float base_disk_gb = 10.0
    RuntimeAttr runtime_default = object {
                                      mem_gb: 16,
                                      disk_gb: ceil(base_disk_gb + (relatedness_size) * 5.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: somalier_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command {
        bcftools view -R ~{bed_file} ~{vcf_uri} -o ~{basename(vcf_uri, '.vcf.gz')+".somalier.subset.vcf.gz"}
        bcftools index -t ~{basename(vcf_uri, '.vcf.gz')+".somalier.subset.vcf.gz"}
    }

    output {
        File subset_vcf = basename(vcf_uri, '.vcf.gz')+".somalier.subset.vcf.gz"
        File subset_vcf_idx = basename(vcf_uri, '.vcf.gz')+".somalier.subset.vcf.gz.tbi"
    }
}

task relatedness_subset {
    input {
        File sites_uri
        File ref_fasta
        File vcf_uri
        File ped_uri
        File ancestry_labels_1kg
        File somalier_1kg_tar
        File sample_file
        String cohort_prefix
        String somalier_docker
        Float relatedness_cutoff
        Boolean infer_ped
        Boolean unknown_flag
        RuntimeAttr? runtime_attr_override
    }

    Float relatedness_size = size(vcf_uri, "GB") 
    Float base_disk_gb = 10.0
    RuntimeAttr runtime_default = object {
                                      mem_gb: 16,
                                      disk_gb: ceil(base_disk_gb + (relatedness_size) * 5.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: somalier_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String infer_string = if infer_ped then "--infer" else ""
    String unknown_flag_str = if unknown_flag then "-u" else ""
    String new_cohort_prefix = basename(sample_file, '.txt')

    command {
        set -euo pipefail

        bcftools index -t ~{vcf_uri}
        somalier extract -d extracted/ --sites ~{sites_uri} -f ~{ref_fasta} ~{vcf_uri}

        somalier relate --ped ~{ped_uri} ~{infer_string} ~{unknown_flag_str} -o ~{new_cohort_prefix} extracted/*.somalier

        tar -xf ~{somalier_1kg_tar}
        somalier ancestry -o ~{new_cohort_prefix} --labels ~{ancestry_labels_1kg} 1kg-somalier/*.somalier ++ extracted/*.somalier
    }

    output {
        File out_samples = new_cohort_prefix + ".samples.tsv" # creates a .ped like vep_annotated_final_vcf with extra QC columns
        File out_pairs = new_cohort_prefix + ".pairs.tsv" # shows IBS for all possible sample pairs
        File out_groups = new_cohort_prefix + ".groups.tsv" # shows pairs of samples above a certain relatedness
        File out_html = new_cohort_prefix + ".html" # interactive html
        File ancestry_html = new_cohort_prefix + ".somalier-ancestry.html"
        File ancestry_out = new_cohort_prefix + ".somalier-ancestry.tsv"
    }
}