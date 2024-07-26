
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
        hl.set(mt.vep.transcript_consequences.Consequence)).size()==0)

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

# Phasing
tmp_ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
tmp_ped.to_csv(f"{prefix}.ped", sep='\t', index=False)
pedigree = hl.Pedigree.read(f"{prefix}.ped", delimiter='\t')

tm = hl.trio_matrix(mt, pedigree, complete_trios=False)
phased_tm = hl.experimental.phase_trio_matrix_by_transmission(tm, call_field='GT', phased_call_field='PBT_GT')

# Output 2.5: OMIM Recessive --> CompHets + Homozygous in probands + XLR in males 
gene_phased_tm = phased_tm.explode_rows(phased_tm.vep.transcript_consequences)
gene_phased_tm = filter_mt(gene_phased_tm)

# XLR only
xlr_phased_tm = gene_phased_tm.filter_rows(gene_phased_tm.vep.transcript_consequences.OMIM_inheritance_code.matches('4'))
xlr_phased_tm = xlr_phased_tm.filter_entries((xlr_phased_tm.proband_entry.GT.is_non_ref()) &
                            (~xlr_phased_tm.is_female))

# filter out calls that couldn't be phased or are hom ref in proband
omim_rec = gene_phased_tm.filter_entries(gene_phased_tm.proband_entry.GT.is_non_ref())

potential_comp_hets = (omim_rec.group_rows_by(omim_rec.vep.transcript_consequences.Feature)
    .aggregate_rows(locus_alleles = hl.agg.collect(omim_rec.row_key))
    .aggregate_entries(proband_PBT_GT = hl.agg.collect(omim_rec.proband_entry.PBT_GT))).result()

potential_comp_hets = potential_comp_hets.filter_rows(
    hl.agg.count_where(hl.set(potential_comp_hets.proband_PBT_GT).size()>1)>0
)

potential_comp_hets = potential_comp_hets.explode_rows(potential_comp_hets.locus_alleles)
potential_comp_hets = potential_comp_hets.key_rows_by(potential_comp_hets.locus_alleles.locus, potential_comp_hets.locus_alleles.alleles, 'Feature')
potential_comp_hets = potential_comp_hets.filter_entries(hl.set(potential_comp_hets.proband_PBT_GT).size()>1)

omim_rec = omim_rec.key_rows_by('locus', 'alleles', omim_rec.vep.transcript_consequences.Feature)
omim_rec = omim_rec.annotate_entries(proband_PBT_GT_set=hl.set(
    potential_comp_hets[omim_rec.row_key, omim_rec.col_key].proband_PBT_GT))

omim_rec = omim_rec.filter_entries(omim_rec.proband_PBT_GT_set.size()>1)
omim_rec_comp_hets = omim_rec.semi_join_rows(potential_comp_hets.rows()).key_rows_by('locus', 'alleles')#.count()

omim_rec_hom_var = gene_phased_tm.filter_entries(gene_phased_tm.proband_entry.GT.is_hom_var())
omim_rec_hom_var = omim_rec_hom_var.filter_rows(hl.agg.count_where(
    hl.is_defined(omim_rec_hom_var.proband_entry.GT))>0).key_rows_by('locus', 'alleles')

xlr_phased_tm = xlr_phased_tm.annotate_rows(variant_type='XLR')
omim_rec_hom_var = omim_rec_hom_var.annotate_rows(variant_type='hom_var')
omim_rec_comp_hets = omim_rec_comp_hets.annotate_rows(variant_type='comphet')

# Output 3: CompHets + Homozygous in probands + XLR in males 
omim_rec_merged = xlr_phased_tm.union_rows(omim_rec_hom_var.drop('Feature', 'proband_PBT_GT_set'))\
.union_rows(omim_rec_comp_hets.drop('Feature', 'proband_PBT_GT_set'))
omim_rec_merged = omim_rec_merged.filter_rows((hl.agg.count_where(hl.is_defined(omim_rec_merged.proband_entry.GT))>0))
# omim_rec_df = omim_rec_merged.entries().to_pandas()
# omim_rec_df['transmission'] = get_transmission(omim_rec_df)

# export OMIM Recessive CompHet + XLR TSV
omim_rec_merged.entries().export(prefix+'_OMIM_recessive_comphet_XLR.tsv.gz', delimiter='\t')
# omim_rec_df.to_csv(prefix+'_OMIM_recessive_comphet_XLR.tsv.gz', sep='\t', index=False)
