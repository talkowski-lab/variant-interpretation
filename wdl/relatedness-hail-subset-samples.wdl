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
        Int chunk_size=100000
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

    call HailPCA {
        input:
        merged_vcf_file=merged_vcf_file,
        cohort_prefix=cohort_prefix,
        bucket_id=bucket_id,
        hail_docker=hail_docker
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
            bucket_id=bucket_id,
            score_table=HailPCA.score_table
        }
    }

    call helpers.mergeResultsPython as mergeRelatednessQC {
            input:
            input_size=size(checkRelatedness.relatedness_qc, 'GB'),
            tsvs=checkRelatedness.relatedness_qc,
            merged_filename=cohort_prefix + '_relatedness_qc.ped',
            hail_docker=hail_docker
    }

    call helpers.mergeResultsPython as mergeKinshipTSV {
            input:
            input_size=size(checkRelatedness.kinship_tsv, 'GB'),
            tsvs=checkRelatedness.kinship_tsv,
            merged_filename=cohort_prefix + '_kinship.tsv',
            hail_docker=hail_docker
    }

    call relatednessHail.plotRelatedness as plotRelatedness {
        input:
        kinship_tsv=mergeKinshipTSV.merged_tsv,
        ped_uri=ped_uri,
        cohort_prefix=cohort_prefix,
        plot_relatedness_script=plot_relatedness_script,
        hail_docker=hail_docker,
        chunk_size=chunk_size
    }

    output {
        File sex_qc_plots = imputeSex.sex_qc_plots
        File ped_sex_qc = imputeSex.ped_sex_qc
        File relatedness_qc = mergeRelatednessQC.merged_tsv
        File kinship_tsv = mergeKinshipTSV.merged_tsv
        File relatedness_plot = plotRelatedness.relatedness_plot
    }
}

task HailPCA {
    input {
        File merged_vcf_file
        String cohort_prefix
        String bucket_id
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(merged_vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
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
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
    cat <<EOF > hail_pca.py
    import hail as hl
    import pandas as pd
    import numpy as np
    import os
    import sys

    merged_vcf_file = sys.argv[1]
    cohort_prefix = sys.argv[2]
    cores = sys.argv[3]
    mem = int(np.floor(float(sys.argv[4])))
    bucket_id = sys.argv[5]

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{mem}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{mem}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    mt = hl.import_vcf(merged_vcf_file, reference_genome='GRCh38', array_elements_required=False, force_bgz=True)
    eigenvalues, score_table, loading_table = hl.hwe_normalized_pca(mt.GT, k=10, compute_loadings=True)
    score_table_file = f"{bucket_id}/hail/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))}/{cohort_prefix}_wes_pca_score_table_som.ht"
    score_table.write(score_table_file, overwrite=True)
    pd.Series([score_table_file]).to_csv('table.txt', header=None, index=False)
    EOF

    python3 hail_pca.py ~{merged_vcf_file} ~{cohort_prefix} ~{cpu_cores} ~{memory} ~{bucket_id}
    >>>

    output {
        String score_table = read_lines('table.txt')[0]
    }
}