version 1.0

# IMPORT
import "Structs.wdl"
import "TinyResolve.wdl"
import "https://raw.githubusercontent.com/broadinstitute/gatk-sv/main/wdl/GetShardInputs.wdl"
#import "https://raw.githubusercontent.com/broadinstitute/gatk-sv/main/wdl/Utils.wdl" as util
import "PEevidence.wdl"

# WORKFLOW DEFINITION
workflow ResolveCTX{
    input{
        File vcf
        String prefix
        String docker_path
        Array[String] samples
        File manta_vcf_tar
        File cytoband
        Array[File] discfile
        Int samples_per_shard = 25
        String sv_pipeline_docker
        String linux_docker
        RuntimeAttr? runtime_attr_resolve
        RuntimeAttr? runtime_attr_untar

        RuntimeAttr? runtime_attr_subset_ctx_vcf
    }

    call CtxVcf2Bed{
        input:
            vcf = vcf,
            prefix = prefix,
            docker_path = docker_path,
            runtime_attr_override = runtime_attr_subset_ctx_vcf
    }

    call TinyResolve.TinyResolve as TinyResolve{
        input:
            samples = samples,
            manta_vcf_tar = manta_vcf_tar,
            cytoband = cytoband,
            discfile = discfile,
            mei_bed = mei_bed,
            samples_per_shard = samples_per_shard,
            sv_pipeline_docker = sv_pipeline_docker,
            linux_docker = linux_docker,
            runtime_attr_resolve = runtime_attr_resolve,
            runtime_attr_untar = runtime_attr_untar
    }

    call mergeMantaVCF{
        input:
            input_vcfs=TinyResolve.tloc_vcf,
            docker = relatedness_docker,
            prefix = prefix,
            runtime_attr_override = runtime_attr_override_merge
    }

#    call ReformatTinyResolve
#    call MergeVCFtinyResolve
#    call PEevidence

    output{
        File ctx_ref_bed = CtxVcf2Bed.ctx_bed
        File merged_manta_raw = mergeMantaVCF.merged_manta
    }
}

# TASK DEFINITIONS
task CtxVcf2Bed {
  input {
    File vcf
    String prefix
    String docker_path
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output{
    File ctx_bed = "~{prefix}.ctx.ref.split.bed"
  }

  command <<<
    set -e
    bcftools view ~{vcf} | grep -E "^#|SVTYPE=CTX" | bgzip -c > ~{prefix}.ctx.vcf.gz
    svtk vcf2bed -i ALL --include-filters ~{prefix}.ctx.vcf.gz - | bgzip -c > ~{prefix}.ctx.bed.gz

    zcat ~{prefix}.ctx.bed.gz | awk '{print $1"\t"$2"\t"$8"\t"$12"\t"$6"\t"$36"\t"$4}'| tail -n+2 > ~{prefix}.ctx.ref.bed

    awk '{split($5,a,",");for(i in a)if(!seen[a[i]]++)print $1"\t"$2"\t"$3"\t"$4"\t"a[i]"\t"$6"\t"$7}' ~{prefix}.ctx.ref.bed > ~{prefix}.ctx.ref.split.bed
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_path
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task mergeMantaVCF{
    input{
        Array [File] input_vcfs
        String docker
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu: 1,
        mem_gb: 12,
        disk_gb: 4,
        boot_disk_gb: 8,
        preemptible: 3,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File merged_manta = "~{prefix}.manta.tloc.vcf.gz"
    }
    command <<<

        bcftools merge -a ~{sep=' ' input_vcfs} | grep -E "^#|SVTYPE=CTX" | bcftools sort | \
        bcftools view - -Oz -o "~{prefix}.manta.tloc.vcf.gz"

    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}