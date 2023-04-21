version 1.0
    
import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks


workflow annotateVcf {

    input {
        File vcf_file
        File fam_ids
        String prefix
        Int records_per_shard
        String variant_interpretation_docker
        String sv_pipeline_updates_docker
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_override_shard_vcf
        RuntimeAttr? runtime_attr_merge
    }

    call MiniTasks.ScatterVcf as SplitVcf {
        input:
            vcf=vcf_file,
            prefix=prefix,
            records_per_shard=records_per_shard,
            sv_pipeline_docker=sv_pipeline_updates_docker,
            runtime_attr_override=runtime_override_shard_vcf
    }

    scatter (shard in SplitVcf.shard_string) {
        call annotate{
            input:
                vcf_file = shard,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_annotate
        }
    }

    call merge {
        input:
            vcf_files = annotate.annotated_vcf,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_merge
    }

    output {
        File annotated_vcf = merge.merged_vcf
    }

}

task annotate{
    input{
        String vcf_file
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File annotated_vcf = "annotated.vcf.gz"
    }

    command {
        set -euo pipefail

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        bcftools annotate -c INFO/COHORT_AC:=INFO/AC,INFO/COHORT_AF:=INFO/AF,INFO/COHORT_AN:=INFO/AN -a ${vcf_file} ${vcf_file} -O z -o annotated.vcf.gz


    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task merge{
    input{
        Array[File] vcf_files
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_files, "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + input_size * 1.5),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File merged_vcf = "merged.vcf.gz"
    }

    command {
        set -euo pipefail

        for vcf in ~{sep=' ' vcf_files}
        do
            export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
            tabix -p vcf $vcf
        done

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        bcftools concat ~{sep=' ' vcf_files} -O z -o merged.vcf.gz
        
    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}
