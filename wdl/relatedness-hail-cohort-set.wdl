version 1.0

import "relatedness-hail.wdl" as relatedness_hail
import "relatedness-hail-subset-samples.wdl" as relatedness_hail_subset_samples
import "mergeVCFs.wdl" as mergeVCFs

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow RelatednessCohortSet {
    input {
        Array[File]? somalier_vcf_files
        File? somalier_vcf_file_
        Array[File] ped_uri
        Array[String] cohort_prefixes
        File sites_uri
        File bed_file
        String merged_filename  # no file extension
        String relatedness_qc_script
        String plot_relatedness_script
        String sex_qc_script
        String sv_base_mini_docker
        String hail_docker
        String bucket_id
        String genome_build
        Int chunk_size=100000
        Int samples_per_chunk=0
        Boolean sort_after_merge=false
        RuntimeAttr? runtime_attr_rename_vcf
        RuntimeAttr? runtime_attr_subset_vcfs
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_impute_sex
        RuntimeAttr? runtime_attr_hail_pca
        RuntimeAttr? runtime_attr_check_relatedness
        RuntimeAttr? runtime_attr_plot_relatedness
        RuntimeAttr? runtime_attr_merge_results
    }

    if (!defined(somalier_vcf_file_)) {
        scatter (pair in zip(cohort_prefixes, select_first([somalier_vcf_files]))) {
            call renameVCFSamples {
                input:
                    vcf_uri=pair.right,
                    cohort_prefix=pair.left,
                    hail_docker=hail_docker,
                    runtime_attr_override=runtime_attr_rename_vcf
            }
        }
        call mergeVCFs.mergeVCFSamples as mergeVCFs {
            input:
            vcf_files=renameVCFSamples.renamed_vcf_file,
            merged_filename=merged_filename,
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_merge_vcfs
        }
    }

    File merged_vcf_file = select_first([somalier_vcf_file_, mergeVCFs.merged_vcf_file])

    call mergePeds {
        input:
            ped_uris=ped_uri,
            cohort_prefixes=cohort_prefixes,
            merged_filename=merged_filename,
            hail_docker=hail_docker,
            input_size=size(ped_uri, 'GB')
    }

    if (samples_per_chunk==0) {
        call relatedness_hail.Relatedness as Relatedness {
            input:
            vep_vcf_files=[],
            somalier_vcf_file_=merged_vcf_file,
            ped_uri=mergePeds.merged_ped_file,
            sites_uri=sites_uri,
            bed_file=bed_file,
            cohort_prefix=merged_filename,
            relatedness_qc_script=relatedness_qc_script,
            plot_relatedness_script=plot_relatedness_script,
            sex_qc_script=sex_qc_script,
            sv_base_mini_docker=sv_base_mini_docker,
            hail_docker=hail_docker,
            bucket_id=bucket_id,
            genome_build=genome_build,
            chunk_size=chunk_size,
            sort_after_merge=sort_after_merge,
            runtime_attr_subset_vcfs=runtime_attr_subset_vcfs,
            runtime_attr_merge_vcfs=runtime_attr_merge_vcfs,
            runtime_attr_impute_sex=runtime_attr_impute_sex,
            runtime_attr_check_relatedness=runtime_attr_check_relatedness,
            runtime_attr_plot_relatedness=runtime_attr_plot_relatedness
        }
    }

    if (samples_per_chunk>0) {
        call relatedness_hail_subset_samples.Relatedness as Relatedness_subsetSamples {
            input:
            vep_vcf_files=[],
            somalier_vcf_file_=merged_vcf_file,
            ped_uri=mergePeds.merged_ped_file,
            sites_uri=sites_uri,
            bed_file=bed_file,
            samples_per_chunk=samples_per_chunk,
            cohort_prefix=merged_filename,
            relatedness_qc_script=relatedness_qc_script,
            plot_relatedness_script=plot_relatedness_script,
            sv_base_mini_docker=sv_base_mini_docker,
            sex_qc_script=sex_qc_script,
            hail_docker=hail_docker,
            bucket_id=bucket_id,
            genome_build=genome_build,
            chunk_size=chunk_size,
            sort_after_merge=sort_after_merge,
            runtime_attr_subset_vcfs=runtime_attr_subset_vcfs,
            runtime_attr_merge_vcfs=runtime_attr_merge_vcfs,
            runtime_attr_impute_sex=runtime_attr_impute_sex,
            runtime_attr_hail_pca=runtime_attr_hail_pca,
            runtime_attr_check_relatedness=runtime_attr_check_relatedness,
            runtime_attr_plot_relatedness=runtime_attr_plot_relatedness,
            runtime_attr_merge_results=runtime_attr_merge_results
        }
    }

    output {
        String cohort_prefix = merged_filename
        File somalier_vcf_file = merged_vcf_file
        File sex_qc_plots = select_first([Relatedness.sex_qc_plots, Relatedness_subsetSamples.sex_qc_plots])
        File ped_sex_qc = select_first([Relatedness.ped_sex_qc, Relatedness_subsetSamples.ped_sex_qc])
        File relatedness_qc = select_first([Relatedness.relatedness_qc, Relatedness_subsetSamples.relatedness_qc])
        File kinship_tsv = select_first([Relatedness.kinship_tsv, Relatedness_subsetSamples.kinship_tsv])
        File relatedness_plot = select_first([Relatedness.relatedness_plot, Relatedness_subsetSamples.relatedness_plot])
    }
}

task mergePeds {
     input {
        Array[String] cohort_prefixes
        Array[String] ped_uris
        String hail_docker
        String merged_filename
        Float input_size
        RuntimeAttr? runtime_attr_override
    }

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
        cat <<EOF > merge_peds.py
        import pandas as pd
        import numpy as np
        import sys

        tsvs = pd.read_csv(sys.argv[1], header=None)[0].tolist()
        cohort_prefixes = pd.read_csv(sys.argv[2], header=None)[0].tolist()
        merged_filename = sys.argv[3] + '.ped'

        dfs = []
        tot = len(tsvs)
        for i, (cohort, tsv) in enumerate(zip(cohort_prefixes, tsvs)):
            print(f"Loading pedigree {i+1}/{tot}...")
            df = pd.read_csv(tsv, sep='\t').iloc[:, :6]
            df.columns = ['FamilyID', 'IndividualID', 'FatherID', 'MotherID', 'Sex', 'Affected']
            df[['FamilyID', 'IndividualID', 'FatherID', 'MotherID']] = df[['FamilyID', 'IndividualID', 'FatherID', 'MotherID']].astype(str)
            df[['FamilyID', 'IndividualID']] = f"{cohort}:" + df[['FamilyID', 'IndividualID']]

            df['FatherID'] = df.FatherID.map({samp: f"{cohort}:" + samp for samp in df[df.FatherID!='0'].FatherID.tolist()} | {'0': '0'})
            df['MotherID'] = df.MotherID.map({samp: f"{cohort}:" + samp for samp in df[df.MotherID!='0'].MotherID.tolist()} | {'0': '0'})

            dfs.append(df)
        merged_df = pd.concat(dfs)
        merged_df.to_csv(merged_filename, sep='\t', index=False)
        EOF

        python3 merge_peds.py ~{write_lines(ped_uris)} ~{write_lines(cohort_prefixes)} ~{merged_filename} > stdout
    >>>

    output {
        File merged_ped_file = "~{merged_filename}.ped" 
    }   
}

task renameVCFSamples {
    input {
        File vcf_uri
        String cohort_prefix
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_uri, 'GB')
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
    cat <<EOF > merge_vcfs.py
    import pandas as pd
    import numpy as np
    import hail as hl
    import sys
    import os

    vcf_uri = sys.argv[1]
    cohort_prefix = sys.argv[2]
    cores = sys.argv[3]  # string
    mem = int(np.floor(float(sys.argv[4])))

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{mem}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{mem}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    mt = hl.import_vcf(vcf_uri, reference_genome='GRCh38', force_bgz=True, call_fields=[], array_elements_required=False)
    mt = mt.key_cols_by()
    mt = mt.annotate_cols(s=hl.str(f"{cohort_prefix}:")+mt.s).key_cols_by('s')
    try:
        # for haploid (e.g. chrY)
        mt = mt.annotate_entries(
            GT = hl.if_else(
                    mt.GT.ploidy == 1, 
                    hl.call(mt.GT[0], mt.GT[0]),
                    mt.GT)
        )
    except:
        pass

    header = hl.get_vcf_metadata(vcf_uri) 
    hl.export_vcf(mt, f"{os.path.basename(vcf_uri).split('.vcf')[0]}_new_sample_IDs.vcf.bgz", metadata=header)
    EOF

    python3 merge_vcfs.py ~{vcf_uri} ~{cohort_prefix} ~{cpu_cores} ~{memory}
    >>>

    String file_ext = if sub(basename(vcf_uri), '.vcf.gz', '')!=basename(vcf_uri) then '.vcf.gz' else '.vcf.bgz'

    output {
        File renamed_vcf_file = "~{basename(vcf_uri, file_ext)}_new_sample_IDs.vcf.bgz"
    }
}
