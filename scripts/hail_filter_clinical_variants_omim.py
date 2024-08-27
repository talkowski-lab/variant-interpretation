
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
am_threshold = float(sys.argv[6])
mpc_threshold = float(sys.argv[7])
gnomad_rec_threshold = float(sys.argv[8])
gnomad_dom_threshold = float(sys.argv[9])
loeuf_v2_threshold = float(sys.argv[10])
loeuf_v4_threshold = float(sys.argv[11])
build = sys.argv[12]

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

    # filter only canonical transcript or MANE PLUS CLINICAL
    mt = mt.filter_rows((mt.vep.transcript_consequences.CANONICAL=='YES') | 
                        (mt.vep.transcript_consequences.MANE_PLUS_CLINICAL!=''))

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

def get_transmission(phased_tm):
    phased_tm = phased_tm.annotate_entries(transmission=hl.if_else(phased_tm.proband_entry.PBT_GT==hl.parse_call('0|0'), 'uninherited',
            hl.if_else(phased_tm.proband_entry.PBT_GT==hl.parse_call('0|1'), 'inherited_from_mother',
                        hl.if_else(phased_tm.proband_entry.PBT_GT==hl.parse_call('1|0'), 'inherited_from_father',
                                hl.or_missing(phased_tm.proband_entry.PBT_GT==hl.parse_call('1|1'), 'inherited_from_both'))))
    )
    return phased_tm

mt = hl.import_vcf(vcf_file, reference_genome=build, force_bgz=True, call_fields=[], array_elements_required=False)

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

gnomad_fields = [x for x in list(mt.vep.transcript_consequences[0]) if 'gnomAD' in x]
mt = mt.annotate_rows(all_csqs=hl.set(hl.flatmap(lambda x: x, mt.vep.transcript_consequences.Consequence)),  
                             gnomad_popmax_af=hl.max([hl.or_missing(hl.array(hl.set(mt.vep.transcript_consequences[gnomad_field]))[0]!='',
                                    hl.float(hl.array(hl.set(mt.vep.transcript_consequences[gnomad_field]))[0])) 
                             for gnomad_field in gnomad_fields]))

# Phasing
tmp_ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
cropped_ped_uri = f"{os.path.basename(ped_uri).split('.ped')[0]}_crop.ped"
tmp_ped.to_csv(cropped_ped_uri, sep='\t', index=False)
pedigree = hl.Pedigree.read(cropped_ped_uri, delimiter='\t')

tm = hl.trio_matrix(mt, pedigree, complete_trios=False)
phased_tm = hl.experimental.phase_trio_matrix_by_transmission(tm, call_field='GT', phased_call_field='PBT_GT')

# Mendel errors
all_errors, per_fam, per_sample, per_variant = hl.mendel_errors(mt['GT'], pedigree)
all_errors_mt = all_errors.key_by().to_matrix_table(row_key=['locus','alleles'], col_key=['s'], row_fields=['fam_id'])
phased_tm = phased_tm.annotate_entries(mendel_code=all_errors_mt[phased_tm.row_key, phased_tm.col_key].mendel_code)

gene_phased_tm = phased_tm.explode_rows(phased_tm.vep.transcript_consequences)
gene_phased_tm = filter_mt(gene_phased_tm)

# Output 2: OMIM Recessive
# OMIM recessive only
omim_rec_gene_phased_tm = gene_phased_tm.filter_rows(
    (gene_phased_tm.vep.transcript_consequences.OMIM_inheritance_code.matches('2')) |    # OMIM recessive
    (gene_phased_tm.vep.transcript_consequences.OMIM_inheritance_code.matches('4')) |  # XLR
    ((hl.is_missing(gene_phased_tm.vep.transcript_consequences.OMIM_inheritance_code)) &  # not OMIM recessive with gnomAD AF and MPC filters
                            ((gene_phased_tm.gnomad_popmax_af<=gnomad_rec_threshold) | (hl.is_missing(gene_phased_tm.gnomad_popmax_af))) &
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
    ((gene_phased_tm.gnomad_popmax_af<=gnomad_dom_threshold) | (hl.is_missing(gene_phased_tm.gnomad_popmax_af))) & 
        ((gene_phased_tm.vep.transcript_consequences.OMIM_inheritance_code.matches('1')) |   # OMIM dominant with gnomAD AF filter
        ((hl.is_missing(gene_phased_tm.vep.transcript_consequences.OMIM_inheritance_code)) &  # not OMIM dominant with LOEUF v2/v4 and MPC filters
            (((gene_phased_tm.info.MPC>=mpc_threshold) | (hl.is_missing(gene_phased_tm.info.MPC))) |
            ((hl.if_else(gene_phased_tm.vep.transcript_consequences.am_pathogenicity=='', 1, 
            hl.float(gene_phased_tm.vep.transcript_consequences.am_pathogenicity))>=am_threshold) &
                ((hl.if_else(gene_phased_tm.vep.transcript_consequences.LOEUF_v2=='', 0, 
                    hl.float(gene_phased_tm.vep.transcript_consequences.LOEUF_v2))<=loeuf_v2_threshold) | 
                (hl.if_else(gene_phased_tm.vep.transcript_consequences.LOEUF_v4=='', 0, 
                    hl.float(gene_phased_tm.vep.transcript_consequences.LOEUF_v4))<=loeuf_v4_threshold)
                )
            )
            )
        )
        )
    )

omim_dom = omim_dom.filter_entries((omim_dom.proband_entry.GT.is_non_ref()) | 
                                   (omim_dom.mother_entry.GT.is_non_ref()) |
                                   (omim_dom.father_entry.GT.is_non_ref()))
omim_dom = omim_dom.filter_rows((hl.agg.count_where(hl.is_defined(omim_dom.proband_entry.GT))>0))
omim_dom = omim_dom.annotate_rows(variant_category='OMIM_dominant')

omim_dom = get_transmission(omim_dom)

# export OMIM Recessive VCF
hl.export_vcf(omim_rec_mt, prefix+'_OMIM_recessive.vcf.bgz', metadata=header)

# export OMIM Dominant TSV
omim_dom.entries().flatten().export(prefix+'_OMIM_dominant.tsv.gz', delimiter='\t')
