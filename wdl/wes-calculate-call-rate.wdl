version 1.0
    
import "scatterVCF.wdl" as scatterVCF
import "mergeSplitVCF.wdl" as mergeSplitVCF
import "mergeVCFs.wdl" as mergeVCFs

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow calculateCallRate {

    input {
        File vcf_file
        File unpadded_intervals_file
        String cohort_prefix
        String bucket_id
        String hail_docker
        String genome_build='GRCh38'
        RuntimeAttr? runtime_attr_infer_platform
    }

    call calculateCallRateMT {
        input:
            vcf_file=vcf_file,
            unpadded_intervals_file=unpadded_intervals_file,
            hail_docker=hail_docker,
            genome_build=genome_build,
            cohort_prefix=cohort_prefix,
            bucket_id=bucket_id,
            runtime_attr_override=runtime_attr_infer_platform
    }

    output {
        String call_rate_mt = calculateCallRateMT.call_rate_mt
    }
}   

task calculateCallRateMT {
    input {
        File vcf_file
        File unpadded_intervals_file
        String cohort_prefix
        String hail_docker
        String genome_build
        String bucket_id
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
    cat <<EOF > call_rate.py
    import hail as hl
    import os
    import sys
    import datetime

    vcf_file = sys.argv[1]
    unpadded_intervals_file = sys.argv[2]
    genome_build = sys.argv[3]
    cohort_prefix = sys.argv[4]
    bucket_id = sys.argv[5]

    hl.init()
    intervals = hl.import_locus_intervals(unpadded_intervals_file, genome_build)
    mt = hl.import_vcf(cohort_df.loc[cohort, 'file'], array_elements_required=False, 
                           reference_genome=genome_build, force_bgz=True, call_fields=[], find_replace=('nul', '.'))
    call_rate_mt = gnomad.sample_qc.platform.compute_callrate_mt(mt, intervals)
    call_rate_mt.write(f"{bucket_id}/hail/call_rate_mt/{cohort_prefix}_call_rate.mt")
    EOF
    python3 call_rate.py ~{vcf_file} ~{unpadded_intervals_file} ~{genome_build} ~{cohort_prefix} ~{bucket_id}
    >>>

    output {
        String call_rate_mt = "~{bucket_id}/hail/call_rate_mt/~{cohort_prefix}_call_rate.mt"
    }
}