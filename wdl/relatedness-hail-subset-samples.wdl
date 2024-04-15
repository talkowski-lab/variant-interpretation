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
        Array[File]? vep_vcf_files
        File? merged_vep_file
        File ped_uri
        File somalier_vcf
        Int samples_per_chunk
        String cohort_prefix
        String relatedness_qc_script
        String plot_relatedness_script
        String sex_qc_script
        String sv_base_mini_docker
        String hail_docker
        String bucket_id
        Boolean sort_after_merge=false
        Int chunk_size=100000
        RuntimeAttr? runtime_attr_subset_vcfs
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_impute_sex
        RuntimeAttr? runtime_attr_hail_pca
        RuntimeAttr? runtime_attr_check_relatedness
        RuntimeAttr? runtime_attr_plot_relatedness
        RuntimeAttr? runtime_attr_merge_results
    }

    if (!defined(merged_vep_file)) {
        
        scatter (vcf_uri in select_first([vep_vcf_files])) {
            String filename = basename(vcf_uri)
            String prefix = if (sub(filename, ".gz", "")!=filename) then basename(filename, ".vcf.gz") else basename(filename, ".vcf.bgz")
            call helpers.subsetVCFs as subsetVCFs {
                input:
                    somalier_vcf=somalier_vcf,
                    vcf_uri=vcf_uri,
                    vcf_idx=vcf_uri+'.tbi',
                    output_name=prefix + '.somalier.subset.vcf.gz',
                    sv_base_mini_docker=sv_base_mini_docker,
                    runtime_attr_override=runtime_attr_subset_vcfs
            }
        }

        call mergeVCFs.mergeVCFs as mergeVCFs {
        input:
            vcf_files=subsetVCFs.subset_vcf,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix,
            sort_after_merge=sort_after_merge,
            runtime_attr_override=runtime_attr_merge_vcfs
        }
    }
    File merged_vcf_file = select_first([merged_vep_file, mergeVCFs.merged_vcf_file])

    call relatednessHail.imputeSex as imputeSex {
        input:
        vcf_uri=merged_vcf_file,
        somalier_vcf=somalier_vcf,
        ped_uri=ped_uri,
        cohort_prefix=cohort_prefix,
        sex_qc_script=sex_qc_script,
        hail_docker=hail_docker,
        runtime_attr_override=runtime_attr_impute_sex
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
        hail_docker=hail_docker,
        runtime_attr_override=runtime_attr_hail_pca
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
            somalier_vcf=somalier_vcf,
            ped_uri=ped_uri,
            cohort_prefix=cohort_prefix,
            relatedness_qc_script=relatedness_qc_script,
            hail_docker=hail_docker,
            bucket_id=bucket_id,
            score_table=HailPCA.score_table,
            runtime_attr_override=runtime_attr_check_relatedness
        }
    }

    call helpers.mergeResultsPython as mergeRelatednessQC {
            input:
            input_size=size(checkRelatedness.relatedness_qc, 'GB'),
            tsvs=checkRelatedness.relatedness_qc,
            merged_filename=cohort_prefix + '_relatedness_qc.ped',
            hail_docker=hail_docker,
            runtime_attr_override=runtime_attr_merge_results
    }

    call helpers.mergeResultsPython as mergeKinshipTSV {
            input:
            input_size=size(checkRelatedness.kinship_tsv, 'GB'),
            tsvs=checkRelatedness.kinship_tsv,
            merged_filename=cohort_prefix + '_kinship.tsv',
            hail_docker=hail_docker,
            runtime_attr_override=runtime_attr_merge_results
    }

    call removeDuplicates {
        input:
        kinship_tsv=mergeKinshipTSV.merged_tsv,
        relatedness_qc=mergeRelatednessQC.merged_tsv,
        hail_docker=hail_docker,
        chunk_size=chunk_size,
        runtime_attr_override=runtime_attr_merge_results
    }

    call relatednessHail.plotRelatedness as plotRelatedness {
        input:
        kinship_tsv=removeDuplicates.kinship_tsv_unique,
        ped_uri=ped_uri,
        cohort_prefix=cohort_prefix,
        plot_relatedness_script=plot_relatedness_script,
        hail_docker=hail_docker,
        chunk_size=chunk_size,
        runtime_attr_override=runtime_attr_plot_relatedness
    }

    output {
        File somalier_vcf_file = merged_vcf_file
        File sex_qc_plots = imputeSex.sex_qc_plots
        File ped_sex_qc = imputeSex.ped_sex_qc
        File relatedness_qc = removeDuplicates.relatedness_qc_unique
        File kinship_tsv = removeDuplicates.kinship_tsv_unique
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
    import datetime

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

    mt = hl.import_vcf(merged_vcf_file, reference_genome='GRCh38', call_fields=[], array_elements_required=False, force_bgz=True)
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

task removeDuplicates {
    input {
        File kinship_tsv
        File relatedness_qc
        String hail_docker
        Int chunk_size
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([kinship_tsv, relatedness_qc], "GB")
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
    cat <<EOF > remove_duplicates.py
    import pandas as pd
    import numpy as np
    import os
    import sys

    kinship_tsv = sys.argv[1]
    relatedness_qc = sys.argv[2]
    chunk_size = int(sys.argv[3])

    chunks = []
    for chunk in pd.read_csv(kinship_tsv, sep='\t', chunksize=chunk_size):
        chunks.append(chunk)
    kinship_df = pd.concat(chunks)
    kinship_df['pair'] = kinship_df[['i','j']].astype(str).agg(lambda lst: ','.join(sorted(lst)), axis=1)
    kinship_df = kinship_df.sort_values('kin', ascending=False).drop_duplicates('pair')

    rel_df = pd.read_csv(relatedness_qc, sep='\t')
    rel_df = rel_df.sort_values(['mother_kin','father_kin'], ascending=False).drop_duplicates('sample_id')

    kinship_df.to_csv(os.path.basename(kinship_tsv), sep='\t', index=False)
    rel_df.to_csv(os.path.basename(relatedness_qc), sep='\t', index=False)
    EOF

    python3 remove_duplicates.py ~{kinship_tsv} ~{relatedness_qc} ~{chunk_size} 
    >>>

    output {
        File kinship_tsv_unique = basename(kinship_tsv)
        File relatedness_qc_unique = basename(relatedness_qc)
    }
}