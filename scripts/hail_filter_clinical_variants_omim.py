
from pyspark.sql import SparkSession
import hail as hl
import numpy as np
import pandas as pd
import sys
import ast
import os

vcf_file = sys.argv[1]
prefix = sys.argv[2]
cores = sys.argv[3]  # string
mem = int(np.floor(float(sys.argv[4])))
ped_uri = sys.argv[5]
ac_threshold = int(sys.argv[6])
gnomad_af_threshold = float(sys.argv[7])
am_threshold = float(sys.argv[8])
mpc_threshold = float(sys.argv[9])
gnomad_rec_threshold = float(sys.argv[10])
gnomad_dom_threshold = float(sys.argv[11])
loeuf_v2_threshold = float(sys.argv[12])
loeuf_v4_threshold = float(sys.argv[13])

hl.init(min_block_size=128, 
        spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                    "spark.driver.cores": "2",
                    "spark.driver.memory": f"{int(np.floor(mem*0.4))}g",
        #             'spark.hadoop.fs.gs.requester.pays.mode': 'CUSTOM',
        #             'spark.hadoop.fs.gs.requester.pays.buckets': 'hail-datasets-us-central1',
        #             'spark.hadoop.fs.gs.requester.pays.project.id': gcp_project,
                    }, 
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

def filter_mt(mt):
    '''
    mt: can be trio matrix (tm) or matrix table (mt) but must be transcript-level, not variant-level
    '''
    # filter by Consequence
    exclude_csqs = ['intergenic_variant', 'upstream_gene_variant', 'downstream_gene_variant',
                    'synonymous_variant', 'coding_sequence_variant', 'sequence_variant']
    mt = mt.filter_rows(hl.set(exclude_csqs).intersection(
        hl.set(mt.vep.transcript_consequences.Consequence)).size()!=hl.set(mt.vep.transcript_consequences.Consequence).size())

    # filter only canonical transcript
    mt = mt.filter_rows(mt.vep.transcript_consequences.CANONICAL=='YES')

    # filter by Impact and splice/noncoding consequence
    splice_vars = ['splice_donor_5th_base_variant', 'splice_region_variant', 'splice_donor_region_variant']
    keep_vars = ['non_coding_transcript_exon_variant']
    mt = mt.filter_rows(
        (hl.set(splice_vars + keep_vars).intersection(
            hl.set(mt.vep.transcript_consequences.Consequence)).size()>0) |
        (hl.array(['HIGH', 'MODERATE']).contains(
        mt.vep.transcript_consequences.IMPACT))
        )
    return mt 

def get_transmission(df):
    '''
    df: trio matrix (tm) phased with PBT_GT converted to Pandas DataFrame
    returns Pandas Series
    ''' 
    return df['proband_entry.PBT_GT'].astype(str).map({
        '0|0': 'uninherited', 
        '0|1': 'inherited_from_mother', 
        '1|0': 'inherited_from_father', 
        '1|1': 'inherited_from_both', 
        'None': 'unknown'
    })

mt = hl.import_vcf(vcf_file, reference_genome='GRCh38', force_bgz=True, call_fields=[], array_elements_required=False)

header = hl.get_vcf_metadata(vcf_file)
csq_columns = header['info']['CSQ']['Description'].split('Format: ')[1].split('|')

# split VEP CSQ string
mt = mt.annotate_rows(vep=mt.info)
transcript_consequences = mt.vep.CSQ.map(lambda x: x.split('\|'))

transcript_consequences_strs = transcript_consequences.map(lambda x: hl.if_else(hl.len(x)>1, hl.struct(**
                                                       {col: x[i] if col!='Consequence' else x[i].split('&')  
                                                        for i, col in enumerate(csq_columns)}), 
                                                        hl.struct(**{col: hl.missing('str') if col!='Consequence' else hl.array([hl.missing('str')])  
                                                        for i, col in enumerate(csq_columns)})))

mt = mt.annotate_rows(vep=mt.vep.annotate(transcript_consequences=transcript_consequences_strs))
mt = mt.annotate_rows(vep=mt.vep.select('transcript_consequences'))

mt = mt.annotate_rows(
    gnomADg_AF=hl.or_missing(hl.array(hl.set(mt.vep.transcript_consequences.gnomADg_AF))[0]!='', 
                                hl.float(hl.array(hl.set(mt.vep.transcript_consequences.gnomADg_AF))[0])),
    gnomADe_AF=hl.or_missing(hl.array(hl.set(mt.vep.transcript_consequences.gnomADe_AF))[0]!='', 
                hl.float(hl.array(hl.set(mt.vep.transcript_consequences.gnomADe_AF))[0])))
mt = mt.annotate_rows(gnomad_af=hl.max([mt.gnomADg_AF, mt.gnomADe_AF]))

# Phasing
tmp_ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
tmp_ped.to_csv(f"{prefix}.ped", sep='\t', index=False)
pedigree = hl.Pedigree.read(f"{prefix}.ped", delimiter='\t')

tm = hl.trio_matrix(mt, pedigree, complete_trios=False)
phased_tm = hl.experimental.phase_trio_matrix_by_transmission(tm, call_field='GT', phased_call_field='PBT_GT')

gene_phased_tm = phased_tm.explode_rows(phased_tm.vep.transcript_consequences)
gene_phased_tm = filter_mt(gene_phased_tm)

# Output 2: OMIM Recessive
# OMIM recessive only
omim_rec_gene_phased_tm = gene_phased_tm.filter_rows(
    (gene_phased_tm.vep.transcript_consequences.OMIM_inheritance_code.matches('2')) |    # OMIM recessive
    (gene_phased_tm.vep.transcript_consequences.OMIM_inheritance_code.matches('4')) |  # XLR
    ((hl.is_missing(gene_phased_tm.vep.transcript_consequences.OMIM_inheritance_code)) &  # not OMIM recessive with gnomAD AF and MPC filters
                            ((gene_phased_tm.gnomad_af<=gnomad_rec_threshold) | (hl.is_missing(gene_phased_tm.gnomad_af))) &
                            ((gene_phased_tm.info.MPC>=mpc_threshold) | (hl.is_missing(gene_phased_tm.info.MPC))) &
                            (hl.if_else(gene_phased_tm.vep.transcript_consequences.am_pathogenicity=='', 1, 
                               hl.float(gene_phased_tm.vep.transcript_consequences.am_pathogenicity))>=am_threshold)
    )
)

omim_rec_gene_phased_tm = (omim_rec_gene_phased_tm.group_rows_by(omim_rec_gene_phased_tm.locus, omim_rec_gene_phased_tm.alleles)
    .aggregate_rows(vep = hl.agg.collect(omim_rec_gene_phased_tm.vep))).result()

fields = list(omim_rec_gene_phased_tm.vep.transcript_consequences[0])
new_csq = omim_rec_gene_phased_tm.vep.transcript_consequences.scan(lambda i, j: 
                                      hl.str('|').join(hl.array([i]))
                                      +','+hl.str('|').join(hl.array([j[col] if col!='Consequence' else 
                                                                  hl.str('&').join(j[col]) 
                                                                  for col in list(fields)])), '')[-1][1:]
omim_rec_gene_phased_tm = omim_rec_gene_phased_tm.annotate_rows(CSQ=new_csq)
omim_rec_mt = mt.semi_join_rows(omim_rec_gene_phased_tm.rows())
omim_rec_mt = omim_rec_mt.annotate_rows(info=omim_rec_mt.info.annotate(CSQ=omim_rec_gene_phased_tm.rows()[omim_rec_mt.row_key].CSQ))

# Output 3: OMIM Dominant
# TODO: temporary hacky, can be removed when VEP rerun with LOEUF HTs fixed
gene_phased_tm = gene_phased_tm.annotate_rows(vep=gene_phased_tm.vep.annotate(transcript_consequences=gene_phased_tm.vep.transcript_consequences.annotate(
    LOEUF_v2=gene_phased_tm.vep.transcript_consequences.LOEUF_v2.replace('null', ''),
    LOEUF_v4=gene_phased_tm.vep.transcript_consequences.LOEUF_v4.replace('null', '')
)))

omim_dom = gene_phased_tm.filter_rows(
        ((gene_phased_tm.gnomad_af<=gnomad_dom_threshold) | (hl.is_missing(gene_phased_tm.gnomad_af))) & 
        ((gene_phased_tm.vep.transcript_consequences.OMIM_inheritance_code.matches('1')) |   # OMIM dominant with gnomAD AF filter
        ((hl.is_missing(gene_phased_tm.vep.transcript_consequences.OMIM_inheritance_code)) &  # not OMIM dominant with LOEUF v2/v4 and MPC filters
            ((gene_phased_tm.info.MPC>=mpc_threshold) | (hl.is_missing(gene_phased_tm.info.MPC))) &
            (hl.if_else(gene_phased_tm.vep.transcript_consequences.am_pathogenicity=='', 1, 
                hl.float(gene_phased_tm.vep.transcript_consequences.am_pathogenicity))>=am_threshold) &
            (hl.if_else(gene_phased_tm.vep.transcript_consequences.LOEUF_v2=='', 0, 
                hl.float(gene_phased_tm.vep.transcript_consequences.LOEUF_v2))<=loeuf_v2_threshold) | 
            (hl.if_else(gene_phased_tm.vep.transcript_consequences.LOEUF_v4=='', 0, 
                hl.float(gene_phased_tm.vep.transcript_consequences.LOEUF_v4))<=loeuf_v4_threshold)
        ))
    )

omim_dom = omim_dom.filter_entries((omim_dom.proband_entry.GT.is_non_ref()) | 
                                   (omim_dom.mother_entry.GT.is_non_ref()) |
                                   (omim_dom.father_entry.GT.is_non_ref()))
omim_dom = omim_dom.filter_rows((hl.agg.count_where(hl.is_defined(omim_dom.proband_entry.GT))>0))
omim_dom = omim_dom.annotate_rows(variant_type='OMIM_dominant')

# export OMIM Recessive VCF
hl.export_vcf(omim_rec_mt, prefix+'_OMIM_recessive.vcf.bgz', metadata=header)

# export OMIM Dominant TSV
omim_dom.entries().flatten().export(prefix+'_OMIM_dominant.tsv.gz', delimiter='\t')
