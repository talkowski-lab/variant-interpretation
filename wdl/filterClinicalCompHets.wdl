version 1.0

import "mergeVCFs.wdl" as mergeVCFs
import "wes-denovo-helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? gpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow filterClinicalCompHets {
    input {
        File? omim_recessive_vcf
        File? sv_filtered_vcf
        File ped_uri
        File omim_uri
        String cohort_prefix
        Array[String] sv_gene_fields = ['PREDICTED_LOF', 'PREDICTED_INTRAGENIC_EXON_DUP']

        String filter_comphets_xlr_hom_var_script

        String hail_docker
        String sv_base_mini_docker

        String genome_build='GRCh38'
        Int families_per_chunk=500
        Boolean mask=false
    }

    # might not work if SV VCF Sample IDs don't match SNV/Indel VCF Sample IDs
    call helpers.splitFamilies as splitFamilies {
        input:
            vcf_file=select_first([omim_recessive_vcf]),
            ped_uri=ped_uri,
            families_per_chunk=families_per_chunk,
            cohort_prefix=cohort_prefix,
            sv_base_mini_docker=sv_base_mini_docker
    } 

    scatter (sample_file in splitFamilies.family_shard_files) {
        if (defined(omim_recessive_vcf)) {
            call helpers.subsetVCFSamplesHail as subsetVCFSamplesSNVIndels {
                input:
                    samples_file=sample_file,
                    vcf_file=select_first([omim_recessive_vcf]),
                    hail_docker=hail_docker,
                    genome_build=genome_build
            }
        }

        if (defined(sv_filtered_vcf)) {
            call helpers.subsetVCFSamplesHail as subsetVCFSamplesSVs {
                input:
                    samples_file=sample_file,
                    vcf_file=select_first([sv_filtered_vcf]),
                    hail_docker=hail_docker,
                    genome_build=genome_build
            }
        }

        call filterCompHetsXLRHomVar {
            input:
                snv_indel_vcf=select_first([subsetVCFSamplesSNVIndels.vcf_subset, 'NA']),
                sv_vcf=select_first([subsetVCFSamplesSVs.vcf_subset, 'NA']),
                ped_uri=ped_uri,
                omim_uri=omim_uri,
                sv_gene_fields=sv_gene_fields,
                filter_comphets_xlr_hom_var_script=filter_comphets_xlr_hom_var_script,
                genome_build=genome_build,
                hail_docker=hail_docker
        }
    }

    call helpers.mergeResultsPython as mergeCompHetsXLRHomVar {
        input:
            tsvs=filterCompHetsXLRHomVar.comphet_xlr_hom_var,
            hail_docker=hail_docker,
            input_size=size(filterCompHetsXLRHomVar.comphet_xlr_hom_var, 'GB'),
            merged_filename=cohort_prefix+'_comp_hets_xlr_hom_var.tsv.gz',
    }

    output {
        File comphet_xlr_hom_var = mergeCompHetsXLRHomVar.merged_tsv
    }
}

task filterCompHetsXLRHomVar {
    input {
        String snv_indel_vcf
        String sv_vcf
        File ped_uri
        String omim_uri

        Array[String] sv_gene_fields
        String genome_build

        String filter_comphets_xlr_hom_var_script
        String hail_docker
        
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size([snv_indel_vcf, sv_vcf], 'GB')
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

    String variant_types_ = if (snv_indel_vcf!='NA') then 'SV_SNV_Indel' else 'SV'
    String variant_types = if (sv_vcf!='NA') then variant_types else 'SNV_Indel'
    String vcf_file = if (variant_types=='SV') then sv_vcf else snv_indel_vcf
    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    String prefix = basename(vcf_file, file_ext) + '_filtered'

    command {
        curl ~{filter_comphets_xlr_hom_var_script} > filter_vcf.py
        python3 filter_vcf.py ~{snv_indel_vcf} ~{sv_vcf} ~{ped_uri} ~{prefix} ~{omim_uri} \
            ~{sep=',' sv_gene_fields} ~{genome_build} ~{cpu_cores} ~{memory} 
    }

    output {
        File comphet_xlr_hom_var = "~{prefix}_{variant_types}_comp_hets_xlr_hom_var.tsv.gz"
    }
}