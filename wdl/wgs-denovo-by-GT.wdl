version 1.0

import "mergeSplitVCF.wdl" as mergeSplitVCF
import "wgs-denovo-step-01.wdl" as step1
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

workflow getDenovoByGT {
    input {
        Array[File] vep_files
        File ped_uri
        File info_header
        Float af_threshold=0.01
        String cohort_prefix
        String vep_hail_docker
        String hail_docker
        String sv_base_mini_docker
        Int shards_per_chunk=10
        Boolean sort_after_merge=false
        Boolean merge_split_vcf=false
        Boolean bad_header=false
        String denovo_snv_indels_gt_script
    }

    String file_ext = if sub(basename(vep_files[0]), '.vcf', '')!=basename(vep_files[0]) then '.vcf' else '.vcf.bgz'
    
    if (merge_split_vcf) {
        call mergeSplitVCF.splitFile as splitVEPFiles {
            input:
                file=write_lines(vep_files),
                shards_per_chunk=shards_per_chunk,
                cohort_prefix=cohort_prefix,
                vep_hail_docker=vep_hail_docker
        }
        # call step1.saveVCFHeader as saveVCFHeaderChunk {
        #     input:
        #         vcf_uri=vep_files[0],
        #         info_header=info_header,
        #         bad_header=bad_header,
        #         sv_base_mini_docker=sv_base_mini_docker
        # }
        scatter (chunk_file in splitVEPFiles.chunks) {        
            call mergeVCFs.mergeVCFs as mergeChunk {
                input:
                    vcf_files=read_lines(chunk_file),
                    sv_base_mini_docker=sv_base_mini_docker,
                    cohort_prefix=basename(chunk_file),
                    sort_after_merge=sort_after_merge
            }
            call denovoByGT as denovoByGTChunk {
                input:
                    vcf_file=mergeChunk.merged_vcf_file,
                    ped_uri=ped_uri,
                    file_ext=file_ext,
                    # header_file=saveVCFHeaderChunk.header_file,
                    af_threshold=af_threshold,
                    vep_hail_docker=vep_hail_docker,
                    denovo_snv_indels_gt_script=denovo_snv_indels_gt_script
            }
        }
        call helpers.mergeResultsPython as mergeChunks {
            input:
                tsvs=denovoByGTChunk.denovo_gt,
                hail_docker=hail_docker,
                input_size=size(denovoByGTChunk.denovo_gt, 'GB'),
                merged_filename=cohort_prefix+'_denovo_GT_AF_filter.tsv.gz'
        }
    }

    if (!merge_split_vcf) {
        # call step1.saveVCFHeader as saveVCFHeader {
        #     input:
        #         vcf_uri=vep_files[0],
        #         info_header=info_header,
        #         bad_header=bad_header,
        #         sv_base_mini_docker=sv_base_mini_docker
        # }
        scatter (vcf_uri in vep_files) {
            call denovoByGT {
                input:
                    vcf_file=vcf_uri,
                    ped_uri=ped_uri,
                    file_ext=file_ext,
                    # header_file=saveVCFHeader.header_file,
                    af_threshold=af_threshold,
                    vep_hail_docker=vep_hail_docker,
                    denovo_snv_indels_gt_script=denovo_snv_indels_gt_script
            }
        }
        call helpers.mergeResultsPython as mergeResults {
            input:
                tsvs=denovoByGT.denovo_gt,
                hail_docker=hail_docker,
                input_size=size(denovoByGT.denovo_gt, 'GB'),
                merged_filename=cohort_prefix+'_denovo_GT_AF_filter.tsv.gz'
        }
    }

    output {
        File denovo_gt = select_first([mergeChunks.merged_tsv, mergeResults.merged_tsv])
    }
}

task denovoByGT {
    input {
        File vcf_file
        File ped_uri
        # File header_file
        String file_ext
        Float af_threshold
        String vep_hail_docker
        String denovo_snv_indels_gt_script
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 10.0
    RuntimeAttr runtime_default = object {
        mem_gb: 8,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: vep_hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command {
        curl ~{denovo_snv_indels_gt_script} > denovo_snv_indels.py
        python3.9 denovo_snv_indels.py ~{vcf_file} ~{ped_uri} ~{af_threshold} \
        ~{cpu_cores} ~{memory} ~{file_ext}
    }

    output {
        File denovo_gt = "~{basename(vcf_file, file_ext)}_denovo_GT_AF_filter.tsv.gz"
    }
}
