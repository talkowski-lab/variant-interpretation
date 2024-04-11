version 1.0

import "relatedness-hail.wdl" as relatedness_hail
import "relatedness-hail-subset-samples.wdl" as relatedness_hail_subset_samples

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
        Array[File] somalier_vcf_file
        Array[File] ped_uri
        File bed_file
        String merged_filename  # no file extension
        String relatedness_qc_script
        String plot_relatedness_script
        String sex_qc_script
        String sv_base_mini_docker
        String hail_docker
        String bucket_id
        Int chunk_size=0
        Int samples_per_chunk=0
        Boolean sort_after_merge=false
        RuntimeAttr? runtime_attr_subset_vcfs
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_impute_sex
        RuntimeAttr? runtime_attr_hail_pca
        RuntimeAttr? runtime_attr_check_relatedness
        RuntimeAttr? runtime_attr_plot_relatedness
        RuntimeAttr? runtime_attr_merge_results
    }

    call mergeVCFs {
        input:
            vcf_files=somalier_vcf_file,
            sv_base_mini_docker=sv_base_mini_docker,
            merged_filename=merged_filename,
            runtime_attr_override=runtime_attr_merge_vcfs
    }
    
    call mergePeds {
        input:
            ped_uris=ped_uri,
            merged_filename=merged_filename,
            hail_docker=hail_docker,
            input_size=size(ped_uri, 'GB')
    }

    if (samples_per_chunk==0) {
        call relatedness_hail.Relatedness as Relatedness {
            input:
            merged_vep_file=mergeVCFs.merged_vcf_file,
            ped_uri=mergePeds.merged_ped_file,
            bed_file=bed_file,
            cohort_prefix=merged_filename,
            relatedness_qc_script=relatedness_qc_script,
            plot_relatedness_script=plot_relatedness_script,
            sex_qc_script=sex_qc_script,
            sv_base_mini_docker=sv_base_mini_docker,
            hail_docker=hail_docker,
            bucket_id=bucket_id,
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
            merged_vep_file=mergeVCFs.merged_vcf_file,
            ped_uri=mergePeds.merged_ped_file,
            bed_file=bed_file,
            samples_per_chunk=samples_per_chunk,
            cohort_prefix=merged_filename,
            relatedness_qc_script=relatedness_qc_script,
            plot_relatedness_script=plot_relatedness_script,
            sv_base_mini_docker=sv_base_mini_docker,
            sex_qc_script=sex_qc_script,
            hail_docker=hail_docker,
            bucket_id=bucket_id,
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
        File somalier_vcf_file = mergeVCFs.merged_vcf_file
        File sex_qc_plots = select_first([Relatedness.sex_qc_plots, Relatedness_subsetSamples.sex_qc_plots])
        File ped_sex_qc = select_first([Relatedness.ped_sex_qc, Relatedness_subsetSamples.ped_sex_qc])
        File relatedness_qc = select_first([Relatedness.relatedness_qc, Relatedness_subsetSamples.relatedness_qc])
        File kinship_tsv = select_first([Relatedness.kinship_tsv, Relatedness_subsetSamples.kinship_tsv])
        File relatedness_plot = select_first([Relatedness.relatedness_plot, Relatedness_subsetSamples.relatedness_plot])
    }
}

task mergeVCFs {
    input {
        Array[File] vcf_files
        String merged_filename
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_files, "GB")
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        VCFS="~{write_lines(vcf_files)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        for vcf in $(cat vcfs_sorted.list);
        do
            tabix $vcf
        done
        bcftools merge --no-version -Oz --file-list vcfs_sorted.list --output ~{merged_filename}_merged.vcf.gz
        
        mkdir -p tmp
        bcftools sort ~{merged_filename}_merged.vcf.gz -Oz --output ~{merged_filename}_sorted.vcf.gz -T tmp/
        bcftools norm ~{merged_filename}_sorted.vcf.gz -m- -Oz --output ~{merged_filename}_merged_sorted.vcf.gz
    >>>

    output {
        File merged_vcf_file = "~{merged_filename}_merged_sorted.vcf.gz"
    }
}

task mergePeds {
     input {
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
        merged_filename = sys.argv[2] + '.ped'

        dfs = []
        tot = len(tsvs)
        for i, tsv in enumerate(tsvs):
            print(f"Loading pedigree {i+1}/{tot}...")
            df = pd.read_csv(tsv, sep='\t').iloc[:, :6]
            df.columns = ['FamilyID', 'IndividualID', 'FatherID', 'MotherID', 'Sex', 'Affected']
            if not df.empty:
                dfs.append(df)
        merged_df = pd.concat(dfs)
        merged_df = merged_df.drop_duplicates('IndividualID')
        merged_df.to_csv(merged_filename, sep='\t', index=False)
        EOF

        python3 merge_peds.py ~{write_lines(ped_uris)} ~{merged_filename} > stdout
    >>>

    output {
        File merged_ped_file = "~{merged_filename}.ped" 
    }   
}

