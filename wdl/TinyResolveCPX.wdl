version 1.0

import "Structs.wdl"
import "GetShardInputs.wdl"
import "Utils.wdl" as util

# this is edited from "gatk-sv/wdl/TinyResolve.wdl"

# Does prelim translocation resolve from raw manta calls
workflow TinyResolve {
  input {
    Array[String] samples         # Sample ID
    File manta_vcf_tar           # tarballed Manta VCFs
    File cytoband
    File cytoband_idx
    Array[File] discfile
    File mei_bed
    Int samples_per_shard = 25
    String sv_pipeline_docker
    String linux_docker
    RuntimeAttr? runtime_attr_resolve
    RuntimeAttr? runtime_attr_untar
  }

  Int num_samples = length(samples)
  Float num_samples_float = num_samples
  Int num_shards = ceil(num_samples_float / samples_per_shard)

  call util.UntarFiles {
    input:
      tar = manta_vcf_tar,
      glob_suffix = ".vcf.gz",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_untar
  }

  scatter (i in range(num_shards)) {
    call GetShardInputs.GetShardInputs as GetShardSamples {
      input:
        items_per_shard = samples_per_shard,
        shard_number = i,
        num_items = num_samples,
        all_items = samples
    }

    call GetShardInputs.GetShardInputs as GetShardDiscfiles {
      input:
        items_per_shard = samples_per_shard,
        shard_number = i,
        num_items = num_samples,
        all_items = discfile
    }

    call GetShardInputs.GetShardInputs as GetShardVcfs {
      input:
        items_per_shard = samples_per_shard,
        shard_number = i,
        num_items = num_samples,
        all_items = UntarFiles.out
    }

    call ResolveManta {
      input:
        raw_vcfs=GetShardVcfs.shard_items, 
        samples=GetShardSamples.shard_items,
        sv_pipeline_docker = sv_pipeline_docker,
        cytoband=cytoband,
        cytoband_idx=cytoband_idx,
        discfile=GetShardDiscfiles.shard_items, 
#        discfile_idx=GetShardDiscfileIndexes.shard_items,
        mei_bed=mei_bed,
        runtime_attr_override=runtime_attr_resolve
    }
  }

  output {
    Array[File] tloc_manta_vcf = flatten(ResolveManta.cpx_vcf)
    Array[File] tloc_manta_vcf_idx = flatten(ResolveManta.cpx_vcf_idx)
    Array[File] tloc_manta_unresolved_vcf = flatten(ResolveManta.cpx_unresolved_vcf)
  }
}


task ResolveManta {
  input {
    Array[File] raw_vcfs ## manta individual vcfs
    Array[String] samples
    File cytoband_idx
    Array[File] discfile ## pe_disc vcfs
#    Array[File] discfile_idx
    File cytoband
    File cytoband_idx
    File mei_bed
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Int num_samples = length(samples)
  Float input_size = size(discfile,"GiB")
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: ceil(10+input_size),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    vcfs=(~{sep=" " raw_vcfs})
    sample_ids=(~{sep=" " samples})
    discfiles=(~{sep=" " discfile})
    for (( i=0; i<~{num_samples}; i++ ));
    do
      vcf=${vcfs[$i]}
      tabix -p vcf $vcf
      sample_id=${sample_ids[$i]}
      pe=${discfiles[$i]}
      tabix -s1 -b2 -e2 $pe
      sample_no=`printf %03d $i`
       
      ## bash /src/variant-interpretation/scripts/mantatloc_check.sh $vcf $pe ${sample_id} ~{mei_bed} ~{cytoband}

      cat <(zcat $std_vcf|egrep ^##) \
        < (echo "##FORMAT=<ID=manta,Number=1,Type=Integer,Description=\"manta genotype\">") \
        < (echo "##INFO=<ID=MEMBERS,Number=.,Type=String,Description=\"IDs of cluster's constituent records.\">") \
        < (echo "##INFO=<ID=EVIDENCE,Number=.,Type=String,Description=\"Classes of random forest support.\">") \
        < (zcat $vcf | egrep -v ^## | awk '{if ($1!~"#")$8=$8";EVIDENCE=PE;MEMBERS="$3; print}' OFS='\t') \
        | bgzip -c > manta.vcf.gz 

      svtk resolve manta.vcf.gz \
      ${sample_id}.manta.complex.vcf \
      --mei-bed $meibed \
      --cytobands $cytobands \
      --discfile $discfile \
      -u ${sample_id}.manta.unresolved.vcf

      bgzip ${sample_id}.manta.complex.vcf
      bgzip ${sample_id}.manta.unresolved.vcf

      mv ${sample_id}.manta.complex.vcf.gz cpx_${sample_no}.${sample_id}.manta.complex.vcf.gz
      tabix -p vcf cpx_${sample_no}.${sample_id}.manta.complex.vcf.gz
    done
  >>>

  output {
    Array[File] tloc_vcf = glob("cpx_*.vcf.gz")
    Array[File] tloc_vcf_idx = glob("cpx_*.vcf.gz.tbi")
    Array[File] tloc_unresolved_vcf = glob("*unresolved.vcf.gz")
  }
  
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
