version 1.0

import "mergeSplitVCF.wdl" as mergeSplitVCF
import "wgs-denovo-step-01.wdl" as step1
import "mergeVCFs.wdl" as mergeVCFs
import "wes-denovo-helpers.wdl" as helpers
import "prioritizeCSQ.wdl" as prioritizeCSQ

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
        File ped_sex_qc
        Float af_threshold=0.01
        String cohort_prefix
        String vep_hail_docker
        String hail_docker
        String sv_base_mini_docker
        Int chunk_size=100000
        Int shards_per_chunk=10
        Boolean sort_after_merge=false
        Boolean merge_split_vcf=false
        String denovo_snv_indels_gt_script
        String prioritize_csq_script
        String sample_column='SAMPLE'
    }

    String file_ext = if sub(basename(vep_files[0]), '.vcf.gz', '')!=basename(vep_files[0]) then '.vcf.gz' else '.vcf.bgz'
    
    if (merge_split_vcf) {
        call mergeSplitVCF.splitFile as splitVEPFiles {
            input:
                file=write_lines(vep_files),
                shards_per_chunk=shards_per_chunk,
                cohort_prefix=cohort_prefix,
                vep_hail_docker=vep_hail_docker
        }

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
                    ped_sex_qc=ped_sex_qc,
                    file_ext=file_ext,
                    af_threshold=af_threshold,
                    vep_hail_docker=vep_hail_docker,
                    denovo_snv_indels_gt_script=denovo_snv_indels_gt_script,
                    prioritize_csq_script=prioritize_csq_script,
                    sample_column=sample_column
            }
        }
        call helpers.mergeResultsPython as mergeChunks {
            input:
                tsvs=denovoByGTChunk.denovo_gt_csq,
                hail_docker=hail_docker,
                input_size=size(denovoByGTChunk.denovo_gt_csq, 'GB'),
                merged_filename=cohort_prefix+'_denovo_GT_AF_filter.tsv.gz'
        }
    }

    if (!merge_split_vcf) {
        scatter (vcf_uri in vep_files) {
            call denovoByGT {
                input:
                    vcf_file=vcf_uri,
                    ped_sex_qc=ped_sex_qc,
                    file_ext=file_ext,
                    af_threshold=af_threshold,
                    vep_hail_docker=vep_hail_docker,
                    denovo_snv_indels_gt_script=denovo_snv_indels_gt_script,
                    prioritize_csq_script=prioritize_csq_script,
                    sample_column=sample_column
            }
        }
        call helpers.mergeResultsPython as mergeResults {
            input:
                tsvs=denovoByGT.denovo_gt_csq,
                hail_docker=hail_docker,
                input_size=size(denovoByGT.denovo_gt_csq, 'GB'),
                merged_filename=cohort_prefix+'_denovo_GT_AF_filter.tsv.gz'
        }
    }

    File denovo_gt_ = select_first([mergeChunks.merged_tsv, mergeResults.merged_tsv])

    call denovoSampleCounts {
        input:
        denovo_gt=denovo_gt_,
        vep_hail_docker=vep_hail_docker,
        sample_column='SAMPLE',
        chunk_size=chunk_size
    }

    output {
        File denovo_gt = denovo_gt_
        File denovo_gt_counts = denovoSampleCounts.denovo_gt_counts
    }
}

task denovoByGT {
    input {
        File vcf_file
        File ped_sex_qc
        String file_ext
        Float af_threshold
        String vep_hail_docker
        String denovo_snv_indels_gt_script
        String prioritize_csq_script    
        String sample_column
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
    
    String denovo_gt = "~{basename(vcf_file, file_ext)}_denovo_GT_AF_filter.tsv.gz"

    command {
        curl ~{denovo_snv_indels_gt_script} > denovo_snv_indels.py
        python3.9 denovo_snv_indels.py ~{vcf_file} ~{ped_sex_qc} ~{af_threshold} \
        ~{cpu_cores} ~{memory} ~{file_ext}

        curl ~{prioritize_csq_script} > prioritize_csq.py
        python3.9 prioritize_csq.py ~{denovo_gt} ~{cpu_cores} ~{memory} ~{sample_column}
    }

    output {
        File denovo_gt_csq = basename(denovo_gt, '.tsv.gz') + '_prioritized_csq.tsv.gz'
    }
}

task denovoSampleCounts {
    input {
        File denovo_gt
        String sample_column
        String vep_hail_docker
        Int chunk_size
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(denovo_gt, "GB")
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

    command <<<
    cat <<EOF > get_counts.py
    import pandas as pd
    import sys
    import os
    import ast

    denovo_gt = sys.argv[1]
    chunk_size = int(sys.argv[2])
    sample_column = sys.argv[3]

    chunks = []
    for chunk in pd.read_csv(denovo_gt, sep='\t', chunksize=chunk_size):
        chunks.append(chunk)
    tm_denovo_df = pd.concat(chunks)
    counts_df = tm_denovo_df.groupby('TYPE')[sample_column].value_counts().reset_index()

    counts_df.to_csv(f"{os.path.basename(denovo_gt).split('.tsv.gz')[0]}_sample_counts.tsv", sep='\t', index=False)
    EOF

    python3.9 get_counts.py ~{denovo_gt} ~{chunk_size} ~{sample_column}
    >>>

    output {
        File denovo_gt_counts = "~{basename(denovo_gt, '.tsv.gz')}_sample_counts.tsv"
    }
}

