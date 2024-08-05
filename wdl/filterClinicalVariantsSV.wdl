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

workflow filterClinicalVariantsSV {
    input {
        File vcf_file
        File clinvar_bed_with_header
        String cohort_prefix
        String genome_build='GRCh38'
        String hail_docker
        String variant_interpretation_docker

        Float bed_overlap_threshold=0.5
    }

    call vcfToBed {
        input:
        vcf_file=vcf_file,
        cohort_prefix=cohort_prefix,
        variant_interpretation_docker=variant_interpretation_docker
    }

    call intersectBed as intersectClinVar {
        input:
        bed_file=vcfToBed.bed_output,
        ref_bed_with_header=clinvar_bed_with_header,
        cohort_prefix=cohort_prefix,
        bed_overlap_threshold=bed_overlap_threshold,
        variant_interpretation_docker=variant_interpretation_docker
    }

    call annotateVCFWithBed as annotateVCFClinVar {
        input:
        vcf_file=vcf_file,
        intersect_bed=intersectClinVar.intersect_bed,
        ref_bed_with_header=clinvar_bed_with_header,
        genome_build=genome_build,
        hail_docker=hail_docker,
        annot_name='ClinVar'
    }

    output {
        File clinvar_vcf = annotateVCFClinVar.filtered_vcf
        File clinvar_vcf_idx = annotateVCFClinVar.filtered_vcf_idx
    }
}

task vcfToBed {
    input{
        File vcf_file
        String cohort_prefix
        String variant_interpretation_docker
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
        docker: variant_interpretation_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -eou pipefail
        svtk vcf2bed --no-samples ~{vcf_file} ~{cohort_prefix}.bed.gz
    >>>

    output {
        File bed_output = "~{cohort_prefix}.bed.gz"
    }
}

task intersectBed {
    input {
        File bed_file
        File ref_bed_with_header
        String cohort_prefix
        String variant_interpretation_docker

        Float bed_overlap_threshold
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([bed_file, ref_bed_with_header], 'GB')
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
        docker: variant_interpretation_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String ref_bed_with_header_str = basename(ref_bed_with_header, '.bed')

    command <<<
        set -eou pipefail
        bedtools intersect -wao -f ~{bed_overlap_threshold} -r -a ~{bed_file} ~{ref_bed_with_header} | bgzip > ~{cohort_prefix}_{ref_bed_with_header_str}.bed.gz
    >>>

    output {
        File intersect_bed = "~{cohort_prefix}_{ref_bed_with_header_str}.bed.gz"
    }
}

task annotateVCFWithBed {
    input {
        File vcf_file
        File intersect_bed
        File ref_bed_with_header
        String annot_name
        String genome_build
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([vcf_file, intersect_bed], 'GB')
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
    cat <<EOF > filter_vcf.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    vcf_file = sys.argv[1]
    intersect_bed = sys.argv[2]
    ref_bed_with_header_uri = sys.argv[3]
    genome_build = sys.argv[4]
    annot_name = sys.argv[5]
    cores = sys.argv[6]
    mem = int(np.floor(float(sys.argv[7])))

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    mt = hl.import_vcf(vcf_file, force_bgz=vcf_file.split('.')[-1] in ['gz', 'bgz'], 
        reference_genome=genome_build, array_elements_required=False, call_fields=[])
    header = hl.get_vcf_metadata(vcf_file)

    overlap_bed = hl.import_table(intersect_bed, force_bgz=True, no_header=True, types={f"f{i}": 'int' for i in [1,2,6,7]})
    overlap_bed = overlap_bed.annotate(f1=overlap_bed.f1 + 1,  # adjust for bed 0-based coordinates
                                    f6=overlap_bed.f6 + 1)

    fields = list(overlap_bed.row)
    overlap_field = fields[-1]

    overlap_bed = overlap_bed.annotate(sv_len=overlap_bed.f2-overlap_bed.f1, 
                        ref_len=overlap_bed.f7-overlap_bed.f6)
    overlap_bed = overlap_bed.annotate(sv_prop=hl.int(overlap_bed[overlap_field]) / overlap_bed.sv_len, 
                        ref_prop=hl.int(overlap_bed[overlap_field]) / overlap_bed.ref_len)

    # use ref_bed_with_header for annotation column/field names
    ref_bed_with_header = hl.import_table(ref_bed_with_header_uri)                                    
    ref_bed_with_header_idx = range(5, len(fields)-1)
    ref_bed_with_header_mapping = {f"f{ref_bed_with_header_idx[i]}": list(ref_bed_with_header.row)[i].lower().replace(' ', '_') 
                                    for i in range(len(ref_bed_with_header_idx))} | {'sv_prop': f"{annot_name}_overlap"}
    overlap_bed = overlap_bed.rename(ref_bed_with_header_mapping)

    overlap_bed = overlap_bed.annotate(locus=hl.locus(overlap_bed.f0, overlap_bed.f1, genome_build))
    overlap_bed = overlap_bed.annotate(alleles=mt.key_rows_by('locus').rows()[overlap_bed.locus].alleles)
    overlap_bed = overlap_bed.key_by('locus','alleles')

    # annotate original VCF
    annot_fields = list(ref_bed_with_header_mapping.values())[3:]
    mt = mt.annotate_rows(info=mt.info.annotate(
        **{field: overlap_bed[mt.row_key][field] for field in annot_fields}))

    for field in annot_fields:
        header['info'][field] = {'Description': '', 'Number': '.', 'Type': 'String'}

    # export filtered and annotated VCF
    hl.export_vcf(mt, os.path.basename(intersect_bed).split('.bed')[0] + '.vcf.bgz', metadata=header, tabix=True)
    EOF

    python3 filter_vcf.py ~{vcf_file} ~{intersect_bed} ~{ref_bed_with_header} ~{genome_build} \
    ~{annot_name} ~{cpu_cores} ~{memory}
    >>>

    output {
        File filtered_vcf = basename(intersect_bed, '.bed.gz') + '.vcf.bgz'
        File filtered_vcf_idx = basename(intersect_bed, '.bed.gz') + '.vcf.bgz.tbi'
    }
}