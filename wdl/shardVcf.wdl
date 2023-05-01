version 1.0
    
import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks


workflow shardVcf {

    input {

        File vcf_file
        Int records_per_shard
        String prefix
        String sv_pipeline_updates_docker
        RuntimeAttr? runtime_override_shard_vcf
    
    }

    call MiniTasks.ScatterVcf as SplitVcf {
        input:
            vcf=subsetVcf.vcf_output,
            prefix=prefix,
            records_per_shard=records_per_shard,
            sv_pipeline_docker=sv_pipeline_updates_docker,
            runtime_attr_override=runtime_override_shard_vcf
    }

    output {
        Array[File] sharded_vcfs = SplitVcf.shards
    }
}