version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow reannotateMPC_AlphaMissense {
    input {
        File vcf_metrics_tsv_final_pu
        File alpha_missense_file

        String loeuf_v2_uri
        String loeuf_v4_uri
        String mpc_ht_uri
        String hail_docker
        String genome_build='GRCh38'
    }

    call reannotateFinalTSV {
        input:
        vcf_metrics_tsv_final_pu=vcf_metrics_tsv_final_pu,
        alpha_missense_file=alpha_missense_file,
        loeuf_v2_uri=loeuf_v2_uri,
        loeuf_v4_uri=loeuf_v4_uri,
        mpc_ht_uri=mpc_ht_uri,
        hail_docker=hail_docker,
        genome_build=genome_build
    }

    output {
        File vcf_metrics_tsv_final_pu_reannot = reannotateFinalTSV.reannotated_tsv
    }
}

task reannotateFinalTSV {
    input {
        File vcf_metrics_tsv_final_pu
        File alpha_missense_file

        String loeuf_v2_uri
        String loeuf_v4_uri
        String mpc_ht_uri
        String hail_docker
        String genome_build

        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_metrics_tsv_final_pu, 'GB')
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
    cat <<EOF > reannotate.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    vcf_metrics_tsv_final_pu = sys.argv[1]
    alpha_missense_file = sys.argv[2]
    mpc_ht_uri = sys.argv[3]
    loeuf_v2_uri = sys.argv[4]
    loeuf_v4_uri = sys.argv[5]
    build = sys.argv[6]
    cores = sys.argv[7]
    mem = int(np.floor(float(sys.argv[8])))

    hl.init(min_block_size=128, 
            local=f"local[*]", 
            spark_conf={
                        "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                        "spark.speculation": 'true'
                        }, 
            tmp_dir="tmp", local_tmpdir="tmp",
                        )

    try:
        ht = hl.import_table(vcf_metrics_tsv_final_pu, force_bgz=vcf_metrics_tsv_final_pu.split('.')[-1]=='gz')
    except:
        ht = hl.import_table(vcf_metrics_tsv_final_pu, force=True)

    if 'locus' in list(ht.row):
        ht = ht.annotate(locus=hl.parse_locus(ht.locus, reference_genome=build))
    else:
        ht = ht.annotate(locus=hl.locus(ht.CHROM, hl.int(ht.POS), reference_genome=build))

    if 'alleles' in list(ht.row):
        ht = ht.annotate(alleles=ht.alleles.replace("\[",'').replace("\]",'').replace("\'",'').replace(' ','').replace('"','').split(','))
    else:
        ht = ht.annotate(alleles=hl.array([ht.REF, ht.ALT]))

    ht = ht.annotate(protein_variant=ht.Protein_position.join(ht.Amino_acids.split('/'))) 

    # AlphaMissense
    am_fields = ['CHROM','POS', 'REF', 'ALT', 'genome', 'uniprot_id', 'transcript_id', 'protein_variant', 'am_pathogenicity', 'am_class']
    am_ht = hl.import_table(alpha_missense_file, comment='#', no_header=True, force_bgz=True)\
        .rename({f"f{i}": am_fields[i] for i in range(len(am_fields))})  # rename fields

    am_ht = am_ht.annotate(locus=hl.locus(am_ht.CHROM, hl.int(am_ht.POS), reference_genome=build),
                alleles=hl.array([am_ht.REF, am_ht.ALT]))

    ht = ht.key_by('locus', 'alleles', 'protein_variant')
    am_ht = am_ht.key_by('locus', 'alleles', 'protein_variant')

    ht = ht.annotate(am_pathogenicity=am_ht[ht.key].am_pathogenicity,
            am_class=am_ht[ht.key].am_class)

    # LOEUF 
    # annotate LOEUF from gnomAD
    existing_loeuf = [field for field in list(ht.row) if 'loeuf' in field.lower()]
    ht = ht.drop(*existing_loeuf)
    loeuf_v2_ht = hl.read_table(loeuf_v2_uri).key_by('transcript')
    loeuf_v4_ht = hl.read_table(loeuf_v4_uri).key_by('transcript')
    ht = ht.key_by('Feature')
    ht = ht.annotate(LOEUF_v2=loeuf_v2_ht[ht.key]['oe_lof_upper'],
                    LOEUF_v2_decile=loeuf_v2_ht[ht.key]['oe_lof_upper_bin'],
                    LOEUF_v4=loeuf_v4_ht[ht.key]['lof.oe_ci.upper'],
                    LOEUF_v4_decile=loeuf_v4_ht[ht.key]['lof.oe_ci.upper_bin_decile'])

    # MPC
    mpc = hl.read_table(mpc_ht_uri).key_by('locus','alleles')
    ht = ht.key_by('locus', 'alleles')
    ht = ht.annotate(MPC=mpc[ht.key].mpc)

    output_filename = os.path.basename(vcf_metrics_tsv_final_pu).split('.tsv')[0] + '.MPC.LOEUF.AlphaMissense.tsv.gz'
    ht.key_by().drop('protein_variant','locus','alleles').export(output_filename)
    EOF

    python3 reannotate.py ~{vcf_metrics_tsv_final_pu} ~{alpha_missense_file} ~{mpc_ht_uri} \
    ~{loeuf_v2_uri} ~{loeuf_v4_uri} ~{genome_build} ~{cpu_cores} ~{memory} 
    >>>

    String file_ext = if sub(basename(vcf_metrics_tsv_final_pu), '.tsv.gz', '')!=basename(vcf_metrics_tsv_final_pu) then '.tsv.gz' else '.tsv'
   
    output {
        File reannotated_tsv = basename(vcf_metrics_tsv_final_pu, file_ext) + '.MPC.LOEUF.AlphaMissense.tsv.gz'
    }
}