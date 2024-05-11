version 1.0

import "Structs.wdl"

workflow Mosaic{
  input{
    String name
    Int? rare_cutoff
    File metrics
    File cutoffs
    File depth_vcf
    File lookup
    File coverage_file
    File coverage_file_idx
    File fam_file
    File median_file
    String sv_pipeline_docker
    String sv_pipeline_rdtest_docker
    File sd_blacklist
    File igl_blacklist

    RuntimeAttr? runtime_attr_mosaic_depth
    RuntimeAttr? runtime_attr_mosaic_potential

  }
  call GetPotential{
   input:
    name=name,
    lookup=lookup,
    metrics=metrics,
    rare_cutoff=rare_cutoff,
    cutoffs=cutoffs,
    depth_vcf=depth_vcf,
    sd_blacklist=sd_blacklist,
    igl_blacklist=igl_blacklist,
    sv_pipeline_docker=sv_pipeline_docker,
    runtime_attr_override=runtime_attr_mosaic_potential
  }
  

  call RdTest{
   input:
    bed=GetPotential.rare,
    coverage_file=coverage_file,
    coverage_file_idx=coverage_file_idx,
    median_file=median_file,
    fam_file=fam_file,
    prefix=name,
    sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
    runtime_attr_override=runtime_attr_mosaic_depth
  }

  output{
    File rare_potential=GetPotential.rare
    File common_potential=GetPotential.common
    File igvplots=RdTest.plots
    File stats=RdTest.stats
  }
}

task GetPotential{
  input{
    String name
    Int? rare_cutoff
    File metrics
    File cutoffs
    File lookup
    File depth_vcf
    File sd_blacklist
    File igl_blacklist
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 8,
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command<<<
    set -euox pipefail
    cut -f 1,2,3,4,36,37 ~{metrics} > phase3-1_7.rd.metrics
    awk '{if ($1~"depth") print}' phase3-1_7.rd.metrics > phase3-1_7.depth.metrics
    # Get the del/dup median separation cutoff from the cutoffs file
    delmed=$(fgrep Depth ~{cutoffs}|fgrep Median|fgrep DEL |cut -f 2)
    dupmed=$(fgrep Depth ~{cutoffs}|fgrep Median|fgrep DUP |cut -f 2)
    delp=$(fgrep Depth ~{cutoffs}|fgrep RD_log_pval|fgrep DEL |cut -f 2)
    dupp=$(fgrep Depth ~{cutoffs}|fgrep RD_log_pval|fgrep DUP |cut -f 2)
    # Find variants that pass p value but not separation
    awk -v delp="$delp" -v delmed="$delmed" '{if ($3=="DEL" && $5<delmed && $6>delp) print}' phase3-1_7.depth.metrics > del.potentialmosaic.txt
    awk -v dupp="$dupp" -v dupmed="$dupmed" '{if ($3=="DUP" && $5<dupmed && $6>dupp) print}' phase3-1_7.depth.metrics> dup.potentialmosaic.txt
    cat del.potentialmosaic.txt dup.potentialmosaic.txt |cut -f1 > potentialmosaic.txt
    tabix -f ~{depth_vcf}
    tabix -H ~{depth_vcf} > head.txt
    zcat ~{depth_vcf} |fgrep -w -f potentialmosaic.txt >body.txt
    cat head.txt body.txt |bgzip -c > test.vcf.gz
    bash /opt/sv-pipeline/04_variant_resolution/scripts/stitch_fragmented_calls.sh -x 1 test.vcf.gz test1.vcf.gz
    bash /opt/sv-pipeline/04_variant_resolution/scripts/stitch_fragmented_calls.sh -x 1 test1.vcf.gz test2.vcf.gz
    bash /opt/sv-pipeline/04_variant_resolution/scripts/stitch_fragmented_calls.sh -x 1 test2.vcf.gz test3.vcf.gz
    svtk vcf2bed test3.vcf.gz ~{name}.potentialmosaic.bed

    ## rare filtering
    while read chr start end id type sample;do
        n=$(zfgrep "$id:" ~{lookup}|cut -f 8)||true
        if [ "$n" -eq "$n" ] ;then
          if [ "$n" -lt ~{rare_cutoff} ]; then
            printf "$chr\t$start\t$end\t$id\t$type\t$sample\n"
          fi
        fi
    done<~{name}.potentialmosaic.bed > ~{name}.potentialmosaic.rare.bed

    echo -e "#chr\tstart\tend\tid\ttype\tsample" > header.bed
    cat header.bed ~{name}.potentialmosaic.bed | bgzip > ~{name}.potentialmosaic.bed.gz
    cat header.bed ~{name}.potentialmosaic.rare.bed | bgzip > ~{name}.potentialmosaic.rare.bed.gz

    zcat < ~{name}.potentialmosaic.rare.bed.gz | awk '{print $1"_"$5"_"$6"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' | sed -e 's/id/name/g' | sed -e 's/type/svtype/g' | sort -k1,1 -k2,2g > ~{name}.depth.ref.bed
    bedtools merge -d 1000 -i ~{name}.depth.ref.bed -delim "," -c 4,5,6 -o collapse > ~{name}.depth.ref.cluster.bed

## refDepth2.R will find sample outliers but will not filter (this is done in postproc)
    Rscript /opt/RdTest/refDepth2.R ~{name}.depth.ref.cluster.bed ~{name}.depth.ref.cluster.filt.bed ~{name}.outliers.txt

## excluding blacklisted calls
    bedtools intersect -v -f 0.5 -wa -wb -a ~{name}.depth.ref.cluster.filt.bed -b ~{sd_blacklist} > ~{name}.depth.ref.cluster.filt.rmSD.bed
    bedtools intersect -v -f 0.5 -wa -wb -a ~{name}.depth.ref.cluster.filt.rmSD.bed -b ~{igl_blacklist} > ~{name}.depth.ref.cluster.filt.rmSD.rmIGL.bed
    cat header.bed ~{name}.depth.ref.cluster.filt.rmSD.rmIGL.bed | sed -e 's/id/name/g' | sed -e 's/svtype/type/g' | bgzip -c > ~{name}.potentialmosaic.rare.bed.gz

## split.sh split clustered SVs
    rm ~{name}.potentialmosaic.rare.bed
    gunzip ~{name}.potentialmosaic.rare.bed.gz
    while read -r line; do \
      chr=$(echo "$line" | cut -f1) \
      start=$(echo "$line" | cut -f2) \
      end=$(echo "$line" | cut -f3) \
      id=$(echo "$line" | cut -f4) \
      type=$(echo "$line" | cut -f5) \
      samples=$(echo "$line" | cut -f6 | tr ',' '\n') \
      
      for sample in $samples; do \
        echo -e "$chr\t$start\t$end\t$id\t$type\t$sample" 
      done
    done < ~{name}.potentialmosaic.rare.bed  > ~{name}.potentialmosaic.rare2.bed
    rm ~{name}.potentialmosaic.rare.bed
    mv ~{name}.potentialmosaic.rare2.bed ~{name}.potentialmosaic.rare.bed
    bgzip ~{name}.potentialmosaic.rare.bed
## fin    

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
  output{
    #File common="~{name}.potentialmosaic.bed.gz"
    File rare = "~{name}.potentialmosaic.rare.bed.gz"
  }
}


# Run rdtest plot
task RdTest {
  input{
    File bed
    String coverage_file
    File coverage_file_idx
    File median_file
    File fam_file
    String prefix
    String sv_pipeline_rdtest_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail

    zcat ~{bed} | tail -n+2 > rdtest.bed
    /opt/RdTest/localize_bincov.sh rdtest.bed ~{coverage_file}
    awk -v OFS="\t" '{print $1,$2,$3,$4,$6,$5}' rdtest.bed > test.bed
    
    ### split test.bed
    while read -r line; do \
      chr=$(echo "$line" | cut -f1) \
      start=$(echo "$line" | cut -f2) \
      end=$(echo "$line" | cut -f3) \
      id=$(echo "$line" | cut -f4) \
      samples=$(echo "$line" | cut -f5 | tr ',' '\n')
      type=$(echo "$line" | cut -f6) \
      
      for sample in $samples; do \
        echo -e "$chr\t$start\t$end\t$id\t$sample\t$type" 
      done
    done < test.bed  > test2.bed
    rm test.bed
    mv test2.bed test.bed
    ### fin

    mkdir plots
    Rscript /opt/RdTest/RdTest.R \
      -b test.bed \
      -n ~{prefix} \
      -c local_coverage.bed.gz \
      -m ~{median_file} \
      -f ~{fam_file} \
      -o plots \
      -p TRUE
    mv plots/~{prefix}.metrics .
    tar -czf mosaic.tar.gz plots/
  >>>

  output {
    File stats = "~{prefix}.metrics"
    File plots= "mosaic.tar.gz"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + (select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + 100) + " HDD"  ### may 11 2024 Increased disk space for cmg and cp cohorts
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_rdtest_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

