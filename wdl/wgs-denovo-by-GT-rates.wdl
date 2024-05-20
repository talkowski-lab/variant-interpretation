version 1.0

import "mergeSplitVCF.wdl" as mergeSplitVCF
import "wgs-denovo-step-01.wdl" as step1
import "mergeVCFs.wdl" as mergeVCFs
import "wes-denovo-helpers.wdl" as helpers
import "prioritizeCSQ.wdl" as prioritizeCSQ
import "annotateMPCandLOEUF.wdl" as annotateMPCandLOEUF

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow getDenovoByGTRates {
    input {
        File denovo_gt
        File ped_sex_qc

        String mpc_dir
        File mpc_chr22_file
        File loeuf_file

        String cohort_prefix
        String hail_docker

        String annotate_mpc_loeuf_script

        Int chunk_size=100000
    }  

    call splitTSV {
        input:
        tsv=denovo_gt,
        chunk_size=chunk_size,
        hail_docker=hail_docker
    }

    scatter (tsv in splitTSV.tsv_shards) {
        call annotateMPCandLOEUF.annotateMPCandLOEUF as annotateMPCandLOEUF {
            input:
                vcf_metrics_tsv=tsv,
                mpc_dir=mpc_dir,
                mpc_chr22_file=mpc_chr22_file,
                loeuf_file=loeuf_file,
                annotate_mpc_loeuf_script=annotate_mpc_loeuf_script,
                hail_docker=hail_docker
        }
    }

    call helpers.mergeResultsPython as mergeAnnotatedTSVs {
        input:
        tsvs=annotateMPCandLOEUF.vcf_metrics_tsv_annot,
        input_size=size(annotateMPCandLOEUF.vcf_metrics_tsv_annot, 'GB'),
        merged_filename=basename(denovo_gt, '.tsv.gz') + '_annot.tsv.gz',
        hail_docker=hail_docker
    }

    call getRates {
        input:
        denovo_gt=mergeAnnotatedTSVs.merged_tsv,
        ped_sex_qc=ped_sex_qc,
        hail_docker=hail_docker,
        cohort_prefix=cohort_prefix,
        chunk_size=chunk_size
    }

    output {
        File denovo_gt_rates = getRates.denovo_gt_rates
    }
}

task splitTSV {
    input {
        File tsv
        Int chunk_size
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(tsv, 'GB')
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
    cat <<EOF > split_tsv.py
    import pandas as pd
    import numpy as np
    import sys
    import os

    uri = sys.argv[1]
    chunk_size = int(sys.argv[2])

    base_filename = os.path.basename(uri).split('.')[0]

    df = pd.concat(pd.read_csv(uri, sep='\t', chunksize=chunk_size))
    for i, sub_df in enumerate(np.array_split(df, chunk_size)):
        sub_df.to_csv(f"{base_filename}_shard_{i}.tsv.gz", sep='\t', index=False)
    EOF

    python3 ~{tsv} ~{chunk_size}
    >>>

    output {
        Array[File] tsv_shards = glob("*_shard_*.tsv.gz")
    }
}

task getRates {
    input {
        File denovo_gt
        File ped_sex_qc
        String hail_docker
        String cohort_prefix
        Int chunk_size
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(denovo_gt, "GB")
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
    
    command <<<
    cat <<EOF > get_counts.py
    import pandas as pd
    import sys
    import os
    import ast

    denovo_gt_uri = sys.argv[1]
    ped_sex_qc = sys.argv[2]
    cohort_prefix = sys.argv[3]
    chunk_size = int(sys.argv[4])

    def get_rates(cohort, n_samples, df, rates_df):
        '''
        df: required columns csq_key, TYPE, LOEUF_tile, MPC
        '''
        rates_df.loc[cohort, 'PTVs'] = df.isPTV.sum() / n_samples
        rates_df.loc[cohort, 'Missense'] = df.isMIS.sum() / n_samples
        rates_df.loc[cohort, 'Missense (MPC at least 2)'] = df[df.MPC>=2].isMIS.sum() / n_samples
        rates_df.loc[cohort, 'Synonymous'] = df.isSYN.sum() / n_samples
        rates_df.loc[cohort, f"PTVs (bottom 3 LOEUF deciles)"] = df[df.LOEUF_tile<=3].isPTV.sum() / n_samples

        for var_type in ['SNV', 'Indel']:
            rates_df.loc[cohort, f"Coding {var_type}s"] = df[df.TYPE==var_type].isCoding.sum() / n_samples
            rates_df.loc[cohort, f"PTVs ({var_type}s)"] = df[df.TYPE==var_type].isPTV.sum() / n_samples
            rates_df.loc[cohort, f"PTVs (bottom 3 LOEUF deciles) ({var_type}s)"] = df[(df.TYPE==var_type)&(df.LOEUF_tile<=3)].isPTV.sum() / n_samples
        return rates_df

    denovo_gt = pd.concat(pd.read_csv(denovo_gt_uri, sep='\t', chunksize=chunk_size))

    ped = pd.read_csv(ped_sex_qc, sep='\t').iloc[:,:6]

    fam_sizes = ped.family_id.value_counts().to_dict()
    fathers = ped[ped.paternal_id!='0'].set_index('sample_id').paternal_id.to_dict()
    mothers = ped[ped.maternal_id!='0'].set_index('sample_id').maternal_id.to_dict()

    def get_sample_role(row):
        if fam_sizes[row.family_id]==1:
            role = 'Singleton'
        elif (row.maternal_id=='0') & (row.paternal_id=='0'):
            if (row.sample_id in fathers.values()):
                role = 'Father'
            elif (row.sample_id in mothers.values()):
                role = 'Mother'
            else:
                role = 'Unknown'
        elif (row.maternal_id=='-9') & (row.paternal_id=='-9'):
            role = 'Unknown'
        else:
            role = 'Proband'
        return role

    ped['role'] = ped.apply(get_sample_role, axis=1)
    cohort_trios = ped.loc[(ped.role=='Proband')&(ped.paternal_id!='0')&(ped.maternal_id!='0'), ['sample_id','paternal_id','maternal_id']]
    n_trios = cohort_trios.shape[0]

    rates_df = pd.DataFrame()
    rates_df = get_rates(cohort_prefix, n_trios, denovo_gt, rates_df)

    rates_df.to_csv(f"{os.path.basename(denovo_gt).split('.tsv.gz')[0]}_rates.tsv", sep='\t')
    EOF

    python3 get_counts.py ~{denovo_gt} ~{ped_sex_qc} ~{cohort_prefix} ~{chunk_size} 
    >>>

    output {
        File denovo_gt_rates = "~{basename(denovo_gt, '.tsv.gz')}_rates.tsv"
    }
}
