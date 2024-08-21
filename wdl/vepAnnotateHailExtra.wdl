version 1.0
    
import "scatterVCF.wdl" as scatterVCF
import "mergeSplitVCF.wdl" as mergeSplitVCF
import "mergeVCFs.wdl" as mergeVCFs
import "annotateNonCoding.wdl" as annotateNonCoding

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow vepAnnotateHailExtra {

    input {
        Array[File] vep_vcf_files

        File revel_file
        File clinvar_vcf_uri
        File omim_uri
        String mpc_ht_uri
        String loeuf_v2_uri
        String loeuf_v4_uri

        String cohort_prefix
        String hail_docker
        String sv_base_mini_docker
        
        String vep_annotate_hail_extra_python_script
        String split_vcf_hail_script

        String genome_build='GRCh38'

        String? spliceAI_snv_uri
        String? spliceAI_indel_uri        
        String? noncoding_bed
        String? gene_list 
        
        RuntimeAttr? runtime_attr_annotate_extra
    }

    scatter (vcf_shard in vep_vcf_files) {

        if (noncoding_bed!='NA') {
            call annotateNonCoding.annotateFromBed as annotateNonCoding {
                input:
                vcf_file=vcf_shard,
                noncoding_bed=select_first([noncoding_bed]),
                hail_docker=hail_docker,
                runtime_attr_override=runtime_attr_annotate_extra
            }
        }
        call annotateExtra {
            input:
                vcf_file=select_first([annotateNonCoding.noncoding_vcf, vcf_shard]),
                vep_annotate_hail_extra_python_script=vep_annotate_hail_extra_python_script,
                loeuf_v2_uri=loeuf_v2_uri,
                loeuf_v4_uri=loeuf_v4_uri,
                revel_file=revel_file,
                revel_file_idx=revel_file+'.tbi',
                clinvar_vcf_uri=clinvar_vcf_uri,
                omim_uri=omim_uri,
                noncoding_bed=select_first([noncoding_bed, 'NA']),
                gene_list=select_first([gene_list, 'NA']),
                mpc_ht_uri=mpc_ht_uri,
                spliceAI_snv_uri=select_first([spliceAI_snv_uri, 'NA']),
                spliceAI_indel_uri=select_first([spliceAI_indel_uri, 'NA']),
                hail_docker=hail_docker,
                genome_build=genome_build,
                runtime_attr_override=runtime_attr_annotate_extra
        }
    }

    output {
        Array[File] annot_vcf_files = annotateExtra.annot_vcf_file
        Array[File] annot_vcf_idx = annotateExtra.annot_vcf_idx
    }
}   

task annotateExtra {
    input {
        File vcf_file
        String loeuf_v2_uri
        String loeuf_v4_uri
        File revel_file
        File revel_file_idx
        File clinvar_vcf_uri
        File omim_uri
        
        String noncoding_bed
        String gene_list

        String mpc_ht_uri
        String spliceAI_snv_uri
        String spliceAI_indel_uri

        String hail_docker
        String genome_build
        String vep_annotate_hail_extra_python_script
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
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String filename = basename(vcf_file)
    String prefix = if (sub(filename, "\\.gz", "")!=filename) then basename(vcf_file, ".vcf.gz") else basename(vcf_file, ".vcf.bgz")
    String vep_annotated_vcf_name = "~{prefix}.annot.vcf.bgz"

    command <<<
        curl ~{vep_annotate_hail_extra_python_script} > annotate.py
        python3 annotate.py -i ~{vcf_file} -o ~{vep_annotated_vcf_name} --cores ~{cpu_cores} --mem ~{memory} \
        --build ~{genome_build} --loeuf-v2 ~{loeuf_v2_uri} --loeuf-v4 ~{loeuf_v4_uri} \
        --mpc ~{mpc_ht_uri} --clinvar ~{clinvar_vcf_uri} --omim ~{omim_uri} \
        --spliceAI-snv ~{spliceAI_snv_uri} --spliceAI-indel ~{spliceAI_indel_uri} \
        --revel ~{revel_file} --genes ~{gene_list} --noncoding-bed ~{noncoding_bed}
        cp $(ls . | grep hail*.log) hail_log.txt
    >>>

    output {
        File annot_vcf_file = vep_annotated_vcf_name
        File annot_vcf_idx = vep_annotated_vcf_name + '.tbi'
        File hail_log = "hail_log.txt"
    }
}