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

workflow filterClinicalCompHets {
    input {
        File? omim_recessive_vcf
        File? clinvar_vcf
        File? sv_filtered_vcf
        File ped_uri
        File omim_uri
        String cohort_prefix
        Array[String] sv_gene_fields = ["PREDICTED_BREAKEND_EXONIC","PREDICTED_COPY_GAIN","PREDICTED_DUP_PARTIAL", 
        "PREDICTED_INTRAGENIC_EXON_DUP","PREDICTED_INTRONIC","PREDICTED_LOF","PREDICTED_MSV_EXON_OVERLAP", 
        "PREDICTED_PARTIAL_EXON_DUP","PREDICTED_PROMOTER","PREDICTED_TSS_DUP","PREDICTED_UTR"]

        String filter_comphets_xlr_hom_var_script

        String hail_docker
        String sv_base_mini_docker

        String genome_build='GRCh38'
        Int families_per_chunk=500
        Boolean mask=false

        RuntimeAttr? runtime_attr_split_families
        RuntimeAttr? runtime_attr_subset_vcfs_snv_indel
        RuntimeAttr? runtime_attr_subset_vcfs_sv
        RuntimeAttr? runtime_attr_filter_comphets
        RuntimeAttr? runtime_attr_merge_results
    }

    if (defined(sv_filtered_vcf)) {
        call addSVSamplesToPed {
            input:
            ped_uri=ped_uri,
            vcf_file=select_first([sv_filtered_vcf]),
            genome_build=genome_build,
            hail_docker=hail_docker,
            runtime_attr_override=runtime_attr_subset_vcfs_sv
        }
    }

    call helpers.splitFamilies as splitFamilies {
        input:
            ped_uri=select_first([addSVSamplesToPed.output_ped, ped_uri]),
            families_per_chunk=families_per_chunk,
            cohort_prefix=cohort_prefix,
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_split_families
    } 

    scatter (sample_file in splitFamilies.family_shard_files) {
        if (defined(omim_recessive_vcf)) {
            call helpers.subsetVCFSamplesHail as subsetVCFSamplesSNVIndels {
                input:
                    samples_file=sample_file,
                    vcf_file=select_first([omim_recessive_vcf]),
                    hail_docker=hail_docker,
                    genome_build=genome_build,
                    runtime_attr_override=runtime_attr_subset_vcfs_snv_indel
            }
        }

        if (defined(clinvar_vcf)) {
            call helpers.subsetVCFSamplesHail as subsetVCFSamplesSNVIndelsClinVar {
                input:
                    samples_file=sample_file,
                    vcf_file=select_first([clinvar_vcf]),
                    hail_docker=hail_docker,
                    genome_build=genome_build,
                    runtime_attr_override=runtime_attr_subset_vcfs_snv_indel
            }
        }

        if (defined(sv_filtered_vcf)) {
            call helpers.subsetVCFSamplesHail as subsetVCFSamplesSVs {
                input:
                    samples_file=sample_file,
                    vcf_file=select_first([sv_filtered_vcf]),
                    hail_docker=hail_docker,
                    genome_build=genome_build,
                    runtime_attr_override=runtime_attr_subset_vcfs_sv
            }
        }

        call filterCompHetsXLRHomVar {
            input:
                snv_indel_vcf=select_first([subsetVCFSamplesSNVIndels.vcf_subset, 'NA']),
                clinvar_vcf=select_first([subsetVCFSamplesSNVIndelsClinVar.vcf_subset, 'NA']),
                sv_vcf=select_first([subsetVCFSamplesSVs.vcf_subset, 'NA']),
                ped_uri=select_first([addSVSamplesToPed.output_ped, ped_uri]),
                omim_uri=omim_uri,
                sv_gene_fields=sv_gene_fields,
                filter_comphets_xlr_hom_var_script=filter_comphets_xlr_hom_var_script,
                genome_build=genome_build,
                hail_docker=hail_docker,
                runtime_attr_override=runtime_attr_filter_comphets
        }
    }

    String variant_types_ = if defined(omim_recessive_vcf) then 'SV_SNV_Indel' else 'SV'
    String variant_types = if defined(sv_filtered_vcf) then variant_types_ else 'SNV_Indel'
    call helpers.mergeResultsPython as mergeCompHetsXLRHomVar {
        input:
            tsvs=filterCompHetsXLRHomVar.comphet_xlr_hom_var_tsv,
            hail_docker=hail_docker,
            input_size=size(filterCompHetsXLRHomVar.comphet_xlr_hom_var_tsv, 'GB'),
            merged_filename="~{cohort_prefix}_~{variant_types}_comp_hets_xlr_hom_var.tsv.gz",
            runtime_attr_override=runtime_attr_merge_results
    }

    output {
        File comphet_xlr_hom_var_tsv = mergeCompHetsXLRHomVar.merged_tsv
    }
}

task addSVSamplesToPed {
    input {
        File ped_uri
        File vcf_file
        String genome_build
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, 'GB')
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
    set -eou pipefail
    cat <<EOF > edit_pedigree.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    vcf_file = sys.argv[1]
    ped_uri = sys.argv[2]
    genome_build = sys.argv[3]
    cores = sys.argv[4]
    mem = int(np.floor(float(sys.argv[5])))

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    mt = hl.import_vcf(vcf_file, force_bgz=vcf_file.split('.')[-1] in ['gz', 'bgz'], 
        reference_genome=genome_build, array_elements_required=False, call_fields=[])
    vcf_samps = mt.s.collect()

    ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
    ped.columns = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
    missing_samps = pd.DataFrame({'family_id': [-9 for _ in range(np.setdiff1d(vcf_samps, ped.sample_id).size)],
            'paternal_id': [0 for _ in range(np.setdiff1d(vcf_samps, ped.sample_id).size)],
            'maternal_id': [0 for _ in range(np.setdiff1d(vcf_samps, ped.sample_id).size)],
            'sample_id': np.setdiff1d(vcf_samps, ped.sample_id)})

    new_ped = pd.concat([ped, missing_samps])
    new_ped = new_ped.replace({np.nan: -9})
    new_ped[['sex','phenotype']] = new_ped[['sex', 'phenotype']].astype(int)
    new_ped.to_csv(os.path.basename(ped_uri).split('.ped')[0] + '_SV_samples.ped',
            sep='\t', index=False)
    EOF
    python3 edit_pedigree.py ~{vcf_file} ~{ped_uri} ~{genome_build} ~{cpu_cores} ~{memory}
    >>>

    output {
        File output_ped = basename(ped_uri, '.ped') + '_SV_samples.ped'
    }
}

task filterCompHetsXLRHomVar {
    input {
        String snv_indel_vcf
        String clinvar_vcf
        String sv_vcf
        File ped_uri
        String omim_uri

        Array[String] sv_gene_fields
        String genome_build

        String filter_comphets_xlr_hom_var_script
        String hail_docker
        
        RuntimeAttr? runtime_attr_override
    }
    String variant_types_ = if (snv_indel_vcf!='NA') then 'SV_SNV_Indel' else 'SV'
    String variant_types = if (sv_vcf!='NA') then variant_types_ else 'SNV_Indel'
    Map[String, Array[String]] vcf_files = {'SV_SNV_Indel': [snv_indel_vcf, clinvar_vcf, sv_vcf], 'SV': [sv_vcf], 'SNV_Indel': [snv_indel_vcf, clinvar_vcf]}

    Float input_size = size(vcf_files[variant_types], 'GB')
    Float base_disk_gb = 10.0
    Float input_disk_scale = 10.0

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

    String vcf_file = if (variant_types=='SV') then sv_vcf else snv_indel_vcf
    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    String prefix = basename(vcf_file, file_ext) + '_filtered'

    command {
        curl ~{filter_comphets_xlr_hom_var_script} > filter_vcf.py
        python3 filter_vcf.py ~{snv_indel_vcf} ~{clinvar_vcf} ~{sv_vcf} ~{ped_uri} ~{prefix} ~{omim_uri} \
            ~{sep=',' sv_gene_fields} ~{genome_build} ~{cpu_cores} ~{memory} 
    }

    output {
        File comphet_xlr_hom_var_tsv = "~{prefix}_~{variant_types}_comp_hets_xlr_hom_var.tsv.gz"
    }
}