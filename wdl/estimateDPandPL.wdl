version 1.0
    
struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow estimateDPandPL {
    input {
        Array[File] vcf_files
        String hail_docker
        String genome_build='GRCh38'
    }

    scatter (vcf_file in vcf_files) {
        call calculateDPandPL {
            input:
            vcf_file=vcf_file,
            hail_docker=hail_docker,
            genome_build=genome_build
        }
    }

    output {
        Array[File] output_vcf_files = calculateDPandPL.output_vcf
        Array[File] output_vcf_idx = calculateDPandPL.output_vcf_idx
    }
}

task calculateDPandPL {
    input {
        File vcf_file
        String hail_docker
        String genome_build
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([vcf_file, vcf_file], "GB") 
    Float base_disk_gb = 10.0

    RuntimeAttr runtime_default = object {
                                      mem_gb: 4,
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String filename = basename(vcf_file)
    String prefix = if (sub(filename, "\\.gz", "")!=filename) then basename(vcf_file, ".vcf.gz") else basename(vcf_file, ".vcf.bgz")
    String output_vcf_name = "~{prefix}.recalc.DP.PL.vcf.bgz"

    command <<<
    cat <<EOF > calculate_dp_pl.py
    from pyspark.sql import SparkSession
    import hail as hl
    import numpy as np
    import pandas as pd
    import sys
    import ast
    import os
    import json
    import argparse

    parser = argparse.ArgumentParser(description='Parse arguments')
    parser.add_argument('-i', dest='vcf_file', help='Input VCF file before annotation')
    parser.add_argument('-o', dest='output_vcf_name', help='Output filename')
    parser.add_argument('--cores', dest='cores', help='CPU cores')
    parser.add_argument('--mem', dest='mem', help='Memory')
    parser.add_argument('--build', dest='build', help='Genome build')

    args = parser.parse_args()

    vcf_file = args.vcf_file
    output_vcf_name = args.output_vcf_name
    cores = args.cores  # string
    mem = int(np.floor(float(args.mem)))
    build = args.build

    hl.init(min_block_size=128, 
            local=f"local[*]", 
            spark_conf={
                        "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                        "spark.speculation": 'true'
                        }, 
            tmp_dir="tmp", local_tmpdir="tmp",
                        )

    #split-multi
    def split_multi_ssc(mt):
        mt = mt.annotate_rows(num_alleles = mt.alleles.size() ) # Add number of alleles at site before split
        # only split variants that aren't already split
        bi = mt.filter_rows(hl.len(mt.alleles) == 2)
        bi = bi.annotate_rows(a_index=1, was_split=False, old_locus=bi.locus, old_alleles=bi.alleles)
        multi = mt.filter_rows(hl.len(mt.alleles) > 2)
        # Now split
        split = hl.split_multi(multi, permit_shuffle=True)
        sm = split.union_rows(bi)
        # sm = hl.split_multi(mt, permit_shuffle=True)
        if 'PL' in list(mt.entry.keys()):
            pl = hl.or_missing(hl.is_defined(sm.PL),
                            (hl.range(0, 3).map(lambda i: hl.min(hl.range(0, hl.len(sm.PL))
            .filter(lambda j: hl.downcode(hl.unphased_diploid_gt_index_call(j), sm.a_index) == hl.unphased_diploid_gt_index_call(i))
            .map(lambda j: sm.PL[j])))))
            sm = sm.annotate_entries(PL = pl)
        split_ds = sm.annotate_entries(GT = hl.downcode(sm.GT, sm.a_index),
                                    AD = hl.or_missing(hl.is_defined(sm.AD), [hl.sum(sm.AD) - sm.AD[sm.a_index], sm.AD[sm.a_index]])
                                    ) 
            #GQ = hl.cond(hl.is_defined(pl[0]) & hl.is_defined(pl[1]) & hl.is_defined(pl[2]), hl.gq_from_pl(pl), sm.GQ) )
        mt = split_ds.drop('old_locus', 'old_alleles')
        return mt

    header = hl.get_vcf_metadata(vcf_file)
    mt = hl.import_vcf(vcf_file, force_bgz=True, array_elements_required=False, call_fields=[], reference_genome=build)
    mt = split_multi_ssc(mt)

    # DP
    mt = mt.annotate_entries(DP=hl.sum(mt.AD))
    header['format']['DP'] = {'Description': 'Approximate read depth (estimated as sum of AD).', 'Number': '1', 'Type': 'Integer'}

    # PL
    mt = mt.annotate_entries(PL=hl.array([0, mt.GQ, hl.max(mt.GQ, 3*(mt.AD[0] - mt.AD[1]))]))
    header['format']['PL'] = {'Description': 'Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification (estimated using GQ).', 'Number': 'G', 'Type': 'Integer'}

    hl.export_vcf(mt, output_vcf_name, metadata=header, tabix=True)
    EOF
    python3 calculate_dp_pl.py -i ~{vcf_file} -o ~{output_vcf_name} --cores ~{cpu_cores} --mem ~{memory} \
        --build ~{genome_build}
    >>>

    output {
        File output_vcf = output_vcf_name
        File output_vcf_idx = output_vcf_name + '.tbi'
    }
}