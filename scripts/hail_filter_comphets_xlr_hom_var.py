from pyspark.sql import SparkSession
import hail as hl
import numpy as np
import pandas as pd
import sys
import ast
import os

from typing import Tuple

import hail.expr.aggregators as agg
from hail.expr import expr_call, expr_float64
from hail.genetics.pedigree import Pedigree
from hail.matrixtable import MatrixTable
from hail.table import Table
from hail.typecheck import numeric, typecheck
from hail.utils.java import Env

snv_indel_vcf = sys.argv[1]
clinvar_vcf = sys.argv[2]
sv_vcf = sys.argv[3]
ped_uri = sys.argv[4]
prefix = sys.argv[5]
omim_uri = sys.argv[6]
sv_gene_fields = sys.argv[7].split(',') 
build = sys.argv[8]
cores = sys.argv[9]  # string
mem = int(np.floor(float(sys.argv[10])))

hl.init(min_block_size=128, 
        local=f"local[*]", 
        spark_conf={
                    "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.speculation": 'true'
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

def load_split_vep_consequences(vcf_uri):
    mt = hl.import_vcf(vcf_uri, reference_genome=build, force_bgz=True, call_fields=[], array_elements_required=False)
    csq_columns = hl.get_vcf_metadata(vcf_uri)['info']['CSQ']['Description'].split('Format: ')[1].split('|')

    mt = mt.annotate_rows(vep=mt.info)
    transcript_consequences = mt.vep.CSQ.map(lambda x: x.split('\|'))

    transcript_consequences_strs = transcript_consequences.map(lambda x: hl.if_else(hl.len(x)>1, hl.struct(**
                                                        {col: x[i] if col!='Consequence' else x[i].split('&')  
                                                            for i, col in enumerate(csq_columns)}), 
                                                            hl.struct(**{col: hl.missing('str') if col!='Consequence' else hl.array([hl.missing('str')])  
                                                            for i, col in enumerate(csq_columns)})))

    mt = mt.annotate_rows(vep=mt.vep.annotate(transcript_consequences=transcript_consequences_strs))
    mt = mt.annotate_rows(vep=mt.vep.select('transcript_consequences'))
    return mt

## STEP 1: Merge SNV/Indel VCF with SV VCF (or just one of them)
# Load SNV/Indel VCF
if snv_indel_vcf!='NA':
    locus_expr = 'locus'
    snv_mt = load_split_vep_consequences(snv_indel_vcf)

    # Load and merge SNV/Indel ClinVar P/LP VCF
    if clinvar_vcf!='NA':
        clinvar_mt = load_split_vep_consequences(clinvar_vcf) 
        snv_mt = snv_mt.union_rows(clinvar_mt).distinct_by_row()

    # filter SNV/Indel MT
    snv_mt = snv_mt.explode_rows(snv_mt.vep.transcript_consequences)
    snv_mt = filter_mt(snv_mt)

    # filter out empty gene fields
    snv_mt = snv_mt.annotate_rows(gene=snv_mt['vep']['transcript_consequences']['SYMBOL'])
    snv_mt = snv_mt.filter_rows(snv_mt.gene!='')

    snv_mt = snv_mt.annotate_rows(variant_type='SNV/Indel')

# Load SV VCF
if sv_vcf!='NA':
    locus_expr = 'locus_interval'
    sv_mt = hl.import_vcf(sv_vcf, reference_genome=build, force_bgz=True, call_fields=[], array_elements_required=False)

    # filter out BNDs
    sv_mt = sv_mt.filter_rows(sv_mt.info.SVTYPE!='BND') 

    # ignore genes for CPX SVs
    sv_mt = sv_mt.annotate_rows(gene=hl.or_missing(sv_mt.info.SVTYPE!='CPX', 
                                                   hl.array(hl.set(hl.flatmap(lambda x: x, [sv_mt.info[field] for field in sv_gene_fields])))))
    sv_mt = sv_mt.annotate_rows(variant_type='SV')

    sv_mt = sv_mt.explode_rows(sv_mt.gene)

    # VEP
    if (snv_indel_vcf!='NA'):
        snv_vep_fields = {field: str(snv_mt.vep.transcript_consequences[field].dtype) for field in list(snv_mt.row.vep.transcript_consequences)}
    else:
        snv_vep_fields = {'OMIM_MIM_number': 'array<str>', 'OMIM_inheritance_code': 'str'}
    sv_mt = sv_mt.annotate_rows(vep=hl.struct(transcript_consequences=
            {field: hl.missing(dtype) for field, dtype in snv_vep_fields.items()}))

    # Annotate OMIM in SVs
    omim = hl.import_table(omim_uri).key_by('approvedGeneSymbol')
    sv_mt = sv_mt.key_rows_by('gene')
    sv_mt = sv_mt.annotate_rows(vep=sv_mt.vep.annotate(
        transcript_consequences=sv_mt.vep.transcript_consequences.annotate(
        OMIM_MIM_number=hl.if_else(hl.is_defined(omim[sv_mt.row_key]), omim[sv_mt.row_key].mimNumber, ''),
        OMIM_inheritance_code=hl.if_else(hl.is_defined(omim[sv_mt.row_key]), omim[sv_mt.row_key].inheritance_code, ''))))
    sv_mt = sv_mt.key_rows_by('locus', 'alleles')

if (snv_indel_vcf!='NA') and (sv_vcf!='NA'):
    sv_info_fields, sv_entry_fields = list(sv_mt.row.info), list(sv_mt.entry)
    snv_info_fields, snv_entry_fields = list(snv_mt.row.info), list(snv_mt.entry)

    sv_missing_entry_fields = {field: str(snv_mt[field].dtype) for field in np.setdiff1d(snv_entry_fields, sv_entry_fields)}
    snv_missing_entry_fields = {field: str(sv_mt[field].dtype) for field in np.setdiff1d(sv_entry_fields, snv_entry_fields)}

    sv_missing_info_fields = {field: str(snv_mt.info[field].dtype) for field in np.setdiff1d(snv_info_fields, sv_info_fields)}
    snv_missing_info_fields = {field: str(sv_mt.info[field].dtype) for field in np.setdiff1d(sv_info_fields, snv_info_fields)}

    sv_mt = sv_mt.annotate_entries(**{field: hl.missing(dtype) for field, dtype in sv_missing_entry_fields.items()})
    snv_mt = snv_mt.annotate_entries(**{field: hl.missing(dtype) for field, dtype in snv_missing_entry_fields.items()})

    sv_mt = sv_mt.select_entries(*sorted(list(sv_mt.entry)))
    snv_mt = snv_mt.select_entries(*sorted(list(snv_mt.entry)))

    sv_mt = sv_mt.annotate_rows(info=sv_mt.info.annotate(**{field: hl.missing(dtype) for field, dtype in sv_missing_info_fields.items()}))
    snv_mt = snv_mt.annotate_rows(info=snv_mt.info.annotate(**{field: hl.missing(dtype) for field, dtype in snv_missing_info_fields.items()}))

    sv_mt = sv_mt.annotate_rows(info=sv_mt.info.select(*sorted(list(sv_mt.info))))
    snv_mt = snv_mt.annotate_rows(info=snv_mt.info.select(*sorted(list(snv_mt.info))))

    sv_mt = sv_mt.key_rows_by().select_rows(*sorted(list(sv_mt.row))).key_rows_by('locus','alleles')
    snv_mt = snv_mt.key_rows_by().select_rows(*sorted(list(snv_mt.row))).key_rows_by('locus','alleles')

    # Subset shared samples 
    sv_samps = sv_mt.s.collect()
    snv_samps = snv_mt.s.collect()
    shared_samps = list(np.intersect1d(sv_samps, snv_samps))

    if len(shared_samps)==0:
        shared_samps = ['']

    def align_mt2_cols_to_mt1(mt1, mt2):
        mt1 = mt1.add_col_index()
        mt2 = mt2.add_col_index()
        new_col_order = mt2.index_cols(mt1.col_key).col_idx.collect()
        return mt2.choose_cols(new_col_order)
    
    sv_mt = sv_mt.filter_cols(hl.array(shared_samps).contains(sv_mt.s))
    snv_mt = align_mt2_cols_to_mt1(sv_mt, snv_mt)

    variant_types = 'SV_SNV_Indel'
    merged_mt = sv_mt.union_rows(snv_mt)
    
    # Change locus to locus_interval to include END for SVs
    merged_mt = merged_mt.annotate_rows(end=hl.if_else(hl.is_defined(merged_mt.info.END2), merged_mt.info.END2, merged_mt.info.END))
    merged_mt = merged_mt.key_rows_by()
    merged_mt = merged_mt.annotate_rows(locus_interval=hl.locus_interval(contig=merged_mt.locus.contig, 
                                                                              start=merged_mt.locus.position,
                                                                              end=merged_mt.end, reference_genome=build))
    locus_expr = 'locus_interval'
    merged_mt = merged_mt.key_rows_by(locus_expr, 'alleles')
    
elif snv_indel_vcf!='NA':
    variant_types = 'SNV_Indel'
    merged_mt = snv_mt

elif sv_vcf!='NA':
    variant_types = 'SV'
    merged_mt = sv_mt

## EDITED HAIL FUNCTIONS
# EDITED: don't check locus struct
@typecheck(dataset=MatrixTable, method=str, tolerate_generic_locus=bool)
def require_biallelic(dataset, method, tolerate_generic_locus: bool = False) -> MatrixTable:
    return dataset._select_rows(
        method,
        hl.case()
        .when(dataset.alleles.length() == 2, dataset._rvrow)
        .or_error(
            f"'{method}' expects biallelic variants ('alleles' field of length 2), found "
            + hl.str(dataset.locus)
            + ", "
            + hl.str(dataset.alleles)
        ),
    )

# EDITED: custom require_biallelic function
@typecheck(call=expr_call, pedigree=Pedigree)
def mendel_errors(call, pedigree) -> Tuple[Table, Table, Table, Table]:
    r"""Find Mendel errors; count per variant, individual and nuclear family.

    .. include:: ../_templates/req_tstring.rst

    .. include:: ../_templates/req_tvariant.rst

    .. include:: ../_templates/req_biallelic.rst

    Examples
    --------

    Find all violations of Mendelian inheritance in each (dad, mom, kid) trio in
    a pedigree and return four tables (all errors, errors by family, errors by
    individual, errors by variant):

    >>> ped = hl.Pedigree.read('data/trios.fam')
    >>> all_errors, per_fam, per_sample, per_variant = hl.mendel_errors(dataset['GT'], ped)

    Export all mendel errors to a text file:

    >>> all_errors.export('output/all_mendel_errors.tsv')

    Annotate columns with the number of Mendel errors:

    >>> annotated_samples = dataset.annotate_cols(mendel=per_sample[dataset.s])

    Annotate rows with the number of Mendel errors:

    >>> annotated_variants = dataset.annotate_rows(mendel=per_variant[dataset.locus, dataset.alleles])

    Notes
    -----

    The example above returns four tables, which contain Mendelian violations
    grouped in various ways. These tables are modeled after the `PLINK mendel
    formats <https://www.cog-genomics.org/plink2/formats#mendel>`_, resembling
    the ``.mendel``, ``.fmendel``, ``.imendel``, and ``.lmendel`` formats,
    respectively.

    **First table:** all Mendel errors. This table contains one row per Mendel
    error, keyed by the variant and proband id.

        - `locus` (:class:`.tlocus`) -- Variant locus, key field.
        - `alleles` (:class:`.tarray` of :py:data:`.tstr`) -- Variant alleles, key field.
        - (column key of `dataset`) (:py:data:`.tstr`) -- Proband ID, key field.
        - `fam_id` (:py:data:`.tstr`) -- Family ID.
        - `mendel_code` (:py:data:`.tint32`) -- Mendel error code, see below.

    **Second table:** errors per nuclear family. This table contains one row
    per nuclear family, keyed by the parents.

        - `pat_id` (:py:data:`.tstr`) -- Paternal ID. (key field)
        - `mat_id` (:py:data:`.tstr`) -- Maternal ID. (key field)
        - `fam_id` (:py:data:`.tstr`) -- Family ID.
        - `children` (:py:data:`.tint32`) -- Number of children in this nuclear family.
        - `errors` (:py:data:`.tint64`) -- Number of Mendel errors in this nuclear family.
        - `snp_errors` (:py:data:`.tint64`) -- Number of Mendel errors at SNPs in this
          nuclear family.

    **Third table:** errors per individual. This table contains one row per
    individual. Each error is counted toward the proband, father, and mother
    according to the `Implicated` in the table below.

        - (column key of `dataset`) (:py:data:`.tstr`) -- Sample ID (key field).
        - `fam_id` (:py:data:`.tstr`) -- Family ID.
        - `errors` (:py:data:`.tint64`) -- Number of Mendel errors involving this
          individual.
        - `snp_errors` (:py:data:`.tint64`) -- Number of Mendel errors involving this
          individual at SNPs.

    **Fourth table:** errors per variant.

        - `locus` (:class:`.tlocus`) -- Variant locus, key field.
        - `alleles` (:class:`.tarray` of :py:data:`.tstr`) -- Variant alleles, key field.
        - `errors` (:py:data:`.tint64`) -- Number of Mendel errors in this variant.

    This method only considers complete trios (two parents and proband with
    defined sex). The code of each Mendel error is determined by the table
    below, extending the
    `Plink classification <https://www.cog-genomics.org/plink2/basic_stats#mendel>`__.

    In the table, the copy state of a locus with respect to a trio is defined
    as follows, where PAR is the `pseudoautosomal region
    <https://en.wikipedia.org/wiki/Pseudoautosomal_region>`__ (PAR) of X and Y
    defined by the reference genome and the autosome is defined by
    :meth:`~.LocusExpression.in_autosome`.

    - Auto -- in autosome or in PAR or female child
    - HemiX -- in non-PAR of X and male child
    - HemiY -- in non-PAR of Y and male child

    `Any` refers to the set \{ HomRef, Het, HomVar, NoCall \} and `~`
    denotes complement in this set.

    +------+---------+---------+--------+----------------------------+
    | Code | Dad     | Mom     | Kid    | Copy State | Implicated    |
    +======+=========+=========+========+============+===============+
    |    1 | HomVar  | HomVar  | Het    | Auto       | Dad, Mom, Kid |
    +------+---------+---------+--------+------------+---------------+
    |    2 | HomRef  | HomRef  | Het    | Auto       | Dad, Mom, Kid |
    +------+---------+---------+--------+------------+---------------+
    |    3 | HomRef  | ~HomRef | HomVar | Auto       | Dad, Kid      |
    +------+---------+---------+--------+------------+---------------+
    |    4 | ~HomRef | HomRef  | HomVar | Auto       | Mom, Kid      |
    +------+---------+---------+--------+------------+---------------+
    |    5 | HomRef  | HomRef  | HomVar | Auto       | Kid           |
    +------+---------+---------+--------+------------+---------------+
    |    6 | HomVar  | ~HomVar | HomRef | Auto       | Dad, Kid      |
    +------+---------+---------+--------+------------+---------------+
    |    7 | ~HomVar | HomVar  | HomRef | Auto       | Mom, Kid      |
    +------+---------+---------+--------+------------+---------------+
    |    8 | HomVar  | HomVar  | HomRef | Auto       | Kid           |
    +------+---------+---------+--------+------------+---------------+
    |    9 | Any     | HomVar  | HomRef | HemiX      | Mom, Kid      |
    +------+---------+---------+--------+------------+---------------+
    |   10 | Any     | HomRef  | HomVar | HemiX      | Mom, Kid      |
    +------+---------+---------+--------+------------+---------------+
    |   11 | HomVar  | Any     | HomRef | HemiY      | Dad, Kid      |
    +------+---------+---------+--------+------------+---------------+
    |   12 | HomRef  | Any     | HomVar | HemiY      | Dad, Kid      |
    +------+---------+---------+--------+------------+---------------+

    See Also
    --------
    :func:`.mendel_error_code`

    Parameters
    ----------
    dataset : :class:`.MatrixTable`
    pedigree : :class:`.Pedigree`

    Returns
    -------
    (:class:`.Table`, :class:`.Table`, :class:`.Table`, :class:`.Table`)
    """
    source = call._indices.source
    if not isinstance(source, MatrixTable):
        raise ValueError(
            "'mendel_errors': expected 'call' to be an expression of 'MatrixTable', found {}".format(
                "expression of '{}'".format(source.__class__) if source is not None else 'scalar expression'
            )
        )

    source = source.select_entries(__GT=call)
    dataset = require_biallelic(source, 'mendel_errors', tolerate_generic_locus=True)
    tm = hl.trio_matrix(dataset, pedigree, complete_trios=True)
    tm = tm.select_entries(
        mendel_code=hl.mendel_error_code(
            tm.locus, tm.is_female, tm.father_entry['__GT'], tm.mother_entry['__GT'], tm.proband_entry['__GT']
        )
    )
    ck_name = next(iter(source.col_key))
    tm = tm.filter_entries(hl.is_defined(tm.mendel_code))
    tm = tm.rename({'id': ck_name})

    entries = tm.entries()

    table1 = entries.select('fam_id', 'mendel_code')

    t2 = tm.annotate_cols(errors=hl.agg.count(), snp_errors=hl.agg.count_where(hl.is_snp(tm.alleles[0], tm.alleles[1])))
    table2 = t2.key_cols_by().cols()
    table2 = table2.select(
        pat_id=table2.father[ck_name],
        mat_id=table2.mother[ck_name],
        fam_id=table2.fam_id,
        errors=table2.errors,
        snp_errors=table2.snp_errors,
    )
    table2 = table2.group_by('pat_id', 'mat_id').aggregate(
        fam_id=hl.agg.take(table2.fam_id, 1)[0],
        children=hl.int32(hl.agg.count()),
        errors=hl.agg.sum(table2.errors),
        snp_errors=hl.agg.sum(table2.snp_errors),
    )
    table2 = table2.annotate(
        errors=hl.or_else(table2.errors, hl.int64(0)), snp_errors=hl.or_else(table2.snp_errors, hl.int64(0))
    )

    # in implicated, idx 0 is dad, idx 1 is mom, idx 2 is child
    implicated = hl.literal(
        [
            [0, 0, 0],  # dummy
            [1, 1, 1],
            [1, 1, 1],
            [1, 0, 1],
            [0, 1, 1],
            [0, 0, 1],
            [1, 0, 1],
            [0, 1, 1],
            [0, 0, 1],
            [0, 1, 1],
            [0, 1, 1],
            [1, 0, 1],
            [1, 0, 1],
        ],
        dtype=hl.tarray(hl.tarray(hl.tint64)),
    )

    table3 = (
        tm.annotate_cols(
            all_errors=hl.or_else(hl.agg.array_sum(implicated[tm.mendel_code]), [0, 0, 0]),
            snp_errors=hl.or_else(
                hl.agg.filter(hl.is_snp(tm.alleles[0], tm.alleles[1]), hl.agg.array_sum(implicated[tm.mendel_code])),
                [0, 0, 0],
            ),
        )
        .key_cols_by()
        .cols()
    )

    table3 = table3.select(
        xs=[
            hl.struct(**{
                ck_name: table3.father[ck_name],
                'fam_id': table3.fam_id,
                'errors': table3.all_errors[0],
                'snp_errors': table3.snp_errors[0],
            }),
            hl.struct(**{
                ck_name: table3.mother[ck_name],
                'fam_id': table3.fam_id,
                'errors': table3.all_errors[1],
                'snp_errors': table3.snp_errors[1],
            }),
            hl.struct(**{
                ck_name: table3.proband[ck_name],
                'fam_id': table3.fam_id,
                'errors': table3.all_errors[2],
                'snp_errors': table3.snp_errors[2],
            }),
        ]
    )
    table3 = table3.explode('xs')
    table3 = table3.select(**table3.xs)
    table3 = (
        table3.group_by(ck_name, 'fam_id')
        .aggregate(errors=hl.agg.sum(table3.errors), snp_errors=hl.agg.sum(table3.snp_errors))
        .key_by(ck_name)
    )

    table4 = tm.select_rows(errors=hl.agg.count_where(hl.is_defined(tm.mendel_code))).rows()

    return table1, table2, table3, table4

## STEP 2: Get CompHets
# Mendel errors
def get_mendel_errors(mt, phased_tm, pedigree):
    all_errors, per_fam, per_sample, per_variant = mendel_errors(mt['GT'], pedigree)  # edited Hail function, see above
    all_errors_mt = all_errors.key_by().to_matrix_table(row_key=[locus_expr,'alleles'], col_key=['s'], col_fields=['fam_id'])
    phased_tm = phased_tm.annotate_entries(mendel_code=all_errors_mt[phased_tm.row_key, phased_tm.col_key].mendel_code)
    return phased_tm

def phase_by_transmission_aggregate_by_gene(tm, mt, pedigree):
    # filter out calls that are hom ref in proband
    tm = tm.filter_entries(tm.proband_entry.GT.is_non_ref())

    phased_tm = hl.experimental.phase_trio_matrix_by_transmission(tm, call_field='GT', phased_call_field='PBT_GT')
    phased_tm = get_mendel_errors(mt, phased_tm, pedigree)
    gene_phased_tm = phased_tm.key_rows_by(locus_expr,'alleles','gene')

    gene_agg_phased_tm = (gene_phased_tm.group_rows_by(gene_phased_tm.gene)
        .aggregate_rows(locus_alleles = hl.agg.collect(gene_phased_tm.row_key),
                       variant_type = hl.agg.collect(gene_phased_tm.variant_type))
        .aggregate_entries(proband_PBT_GT = hl.agg.collect(gene_phased_tm.proband_entry.PBT_GT).filter(hl.is_defined),
                          proband_GT = hl.agg.collect(gene_phased_tm.proband_entry.GT).filter(hl.is_defined))).result()
    return gene_phased_tm, gene_agg_phased_tm

def get_subset_tm(mt, samples, pedigree, keep=True, complete_trios=False):
    subset_mt = mt.filter_cols(hl.array(samples).contains(mt.s), keep=keep)

    # remove variants missing in subset samples
    subset_mt = hl.variant_qc(subset_mt)
    subset_mt = subset_mt.filter_rows(subset_mt.variant_qc.AC[1]>0)
    subset_mt = subset_mt.drop('variant_qc')

    subset_tm = hl.trio_matrix(subset_mt, pedigree, complete_trios=complete_trios)
    return subset_mt, subset_tm

def get_non_trio_comphets(mt):
    non_trio_mt, non_trio_tm = get_subset_tm(mt, non_trio_samples, non_trio_pedigree)
    non_trio_gene_phased_tm, non_trio_gene_agg_phased_tm = phase_by_transmission_aggregate_by_gene(non_trio_tm, non_trio_mt, non_trio_pedigree)

    # different criteria for non-trios
    potential_comp_hets_non_trios = non_trio_gene_agg_phased_tm.filter_rows(
            hl.agg.count_where(non_trio_gene_agg_phased_tm.proband_GT.size()>1)>0
    )
          
    potential_comp_hets_non_trios = potential_comp_hets_non_trios.explode_rows(potential_comp_hets_non_trios.locus_alleles)
    potential_comp_hets_non_trios = potential_comp_hets_non_trios.key_rows_by(potential_comp_hets_non_trios.locus_alleles[locus_expr], potential_comp_hets_non_trios.locus_alleles.alleles, 'gene')

    potential_comp_hets_non_trios = potential_comp_hets_non_trios.filter_entries(potential_comp_hets_non_trios.proband_GT.size()>1)
    
    non_trio_gene_phased_tm = non_trio_gene_phased_tm.key_rows_by(locus_expr, 'alleles', 'gene')
    non_trio_gene_phased_tm = non_trio_gene_phased_tm.annotate_entries(proband_GT=
        potential_comp_hets_non_trios[non_trio_gene_phased_tm.row_key, non_trio_gene_phased_tm.col_key].proband_GT,
                                                                      proband_GT_set=hl.set(
        potential_comp_hets_non_trios[non_trio_gene_phased_tm.row_key, non_trio_gene_phased_tm.col_key].proband_GT),
                                                                      proband_PBT_GT_set=hl.set(
        potential_comp_hets_non_trios[non_trio_gene_phased_tm.row_key, non_trio_gene_phased_tm.col_key].proband_PBT_GT))
    non_trio_gene_phased_tm = non_trio_gene_phased_tm.filter_entries(non_trio_gene_phased_tm.proband_GT.size()>1)  # this actually seems necessary
    gene_phased_tm_comp_hets_non_trios = non_trio_gene_phased_tm.semi_join_rows(potential_comp_hets_non_trios.rows()).key_rows_by(locus_expr, 'alleles')
    return gene_phased_tm_comp_hets_non_trios

def get_trio_comphets(mt):
    trio_mt, trio_tm = get_subset_tm(mt, trio_samples, trio_pedigree, keep=True, complete_trios=True)
    trio_gene_phased_tm, trio_gene_agg_phased_tm = phase_by_transmission_aggregate_by_gene(trio_tm, trio_mt, trio_pedigree)

    # different criteria for trios (requires phasing)
    potential_comp_hets_trios = trio_gene_agg_phased_tm.filter_rows(
        hl.agg.count_where(hl.set(trio_gene_agg_phased_tm.proband_PBT_GT).size()>1)>0
    )
    potential_comp_hets_trios = potential_comp_hets_trios.explode_rows(potential_comp_hets_trios.locus_alleles)
    potential_comp_hets_trios = potential_comp_hets_trios.key_rows_by(potential_comp_hets_trios.locus_alleles[locus_expr], potential_comp_hets_trios.locus_alleles.alleles, 'gene')

    potential_comp_hets_trios = potential_comp_hets_trios.filter_entries(hl.set(potential_comp_hets_trios.proband_PBT_GT).size()>1)

    trio_gene_phased_tm = trio_gene_phased_tm.key_rows_by(locus_expr, 'alleles', 'gene')
    trio_gene_phased_tm = trio_gene_phased_tm.annotate_entries(proband_GT=
        potential_comp_hets_trios[trio_gene_phased_tm.row_key, trio_gene_phased_tm.col_key].proband_GT,
                                                               proband_GT_set=hl.set(
        potential_comp_hets_trios[trio_gene_phased_tm.row_key, trio_gene_phased_tm.col_key].proband_GT),
                                                               proband_PBT_GT_set=hl.set(
        potential_comp_hets_trios[trio_gene_phased_tm.row_key, trio_gene_phased_tm.col_key].proband_PBT_GT))

    trio_gene_phased_tm = trio_gene_phased_tm.filter_entries(trio_gene_phased_tm.proband_PBT_GT_set.size()>1)  # this actually seems necessary
    gene_phased_tm_comp_hets_trios = trio_gene_phased_tm.semi_join_rows(potential_comp_hets_trios.rows()).key_rows_by(locus_expr, 'alleles')
    return gene_phased_tm_comp_hets_trios

def get_transmission(phased_tm_ht):
    phased_tm_ht = phased_tm_ht.annotate(transmission=hl.if_else(phased_tm_ht.proband_entry.PBT_GT==hl.parse_call('0|0'), 'uninherited',
            hl.if_else(phased_tm_ht.proband_entry.PBT_GT==hl.parse_call('0|1'), 'inherited_from_mother',
                        hl.if_else(phased_tm_ht.proband_entry.PBT_GT==hl.parse_call('1|0'), 'inherited_from_father',
                                hl.or_missing(phased_tm_ht.proband_entry.PBT_GT==hl.parse_call('1|1'), 'inherited_from_both'))))
    )
    return phased_tm_ht

# Subset pedigree to samples in VCF, edit parental IDs
vcf_samples = merged_mt.s.collect()
tmp_ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
tmp_ped = tmp_ped[tmp_ped.sample_id.isin(vcf_samples)].copy()
tmp_ped['paternal_id'] = tmp_ped.paternal_id.apply(lambda id: id if id in vcf_samples else '0')
tmp_ped['maternal_id'] = tmp_ped.maternal_id.apply(lambda id: id if id in vcf_samples else '0')
tmp_ped.to_csv(f"{prefix}.ped", sep='\t', index=False)

pedigree = hl.Pedigree.read(f"{prefix}.ped", delimiter='\t')
pedigree = pedigree.filter_to(vcf_samples)

trio_samples = list(np.intersect1d(vcf_samples,
                              list(np.array([[trio.s, trio.pat_id, trio.mat_id] 
                                             for trio in pedigree.complete_trios() if trio.fam_id!='-9']).flatten())))
# Make sure incomplete trios include parents
fathers = pd.Series({trio.s: trio.pat_id for trio in pedigree.trios if trio.pat_id is not None})
mothers = pd.Series({trio.s: trio.mat_id for trio in pedigree.trios if trio.mat_id is not None})
non_trio_samples = list(np.setdiff1d(vcf_samples, trio_samples))
non_trio_samples = list(np.union1d(np.union1d(non_trio_samples,
                              mothers.loc[np.intersect1d(non_trio_samples, mothers.index)].tolist()),
                             fathers.loc[np.intersect1d(non_trio_samples, fathers.index)].tolist()))

trio_pedigree = pedigree.filter_to(trio_samples)
non_trio_pedigree = pedigree.filter_to(non_trio_samples)

## Get CompHets 
# Filter to only in autosomes or PAR
comphet_mt = merged_mt.filter_rows(merged_mt.locus.in_autosome_or_par())
# Filter only OMIM recessive or missing (for SVs)
comphet_mt = comphet_mt.filter_rows((comphet_mt.vep.transcript_consequences.OMIM_inheritance_code.matches('2')) | 
                                    (comphet_mt.vep.transcript_consequences.OMIM_inheritance_code==''))

merged_trio_comphets = get_trio_comphets(comphet_mt)
merged_non_trio_comphets = get_non_trio_comphets(comphet_mt)
merged_trio_comphets = merged_trio_comphets.annotate_cols(trio_status='trio')
merged_non_trio_comphets = merged_non_trio_comphets.annotate_cols(trio_status=
                                                              hl.if_else(merged_non_trio_comphets.fam_id=='-9', 
                                                              'not_in_pedigree', 'non_trio'))
merged_comphets = merged_trio_comphets.entries().union(merged_non_trio_comphets.entries())

# XLR only
merged_tm = hl.trio_matrix(merged_mt, pedigree, complete_trios=False)
gene_phased_tm, gene_agg_phased_tm = phase_by_transmission_aggregate_by_gene(merged_tm, merged_mt, pedigree)
gene_phased_tm = gene_phased_tm.annotate_cols(trio_status=hl.if_else(gene_phased_tm.fam_id=='-9', 'not_in_pedigree', 
                                                   hl.if_else(hl.array(trio_samples).contains(gene_phased_tm.id), 'trio', 'non_trio')))

xlr_phased_tm = gene_phased_tm.filter_rows((gene_phased_tm.vep.transcript_consequences.OMIM_inheritance_code.matches('4')) |  # OMIM XLR
                                           ((gene_phased_tm.locus.in_x_nonpar()) | (gene_phased_tm.locus.in_x_par())))  # on X chromosome
xlr_phased = xlr_phased_tm.filter_entries((xlr_phased_tm.proband_entry.GT.is_non_ref()) &
                            (~xlr_phased_tm.is_female)).key_rows_by(locus_expr, 'alleles').entries()


# HomVar in proband only
phased_hom_var = gene_phased_tm.filter_entries(gene_phased_tm.proband_entry.GT.is_hom_var())
phased_hom_var = phased_hom_var.filter_entries((phased_hom_var.locus.in_x_nonpar()) &
                            (~phased_hom_var.is_female), keep=False)  # filter out non-PAR chrX in males
phased_hom_var = phased_hom_var.filter_rows(hl.agg.count_where(
    hl.is_defined(phased_hom_var.proband_entry.GT))>0).key_rows_by(locus_expr, 'alleles').entries()

xlr_phased = xlr_phased.annotate(variant_category='XLR')
phased_hom_var = phased_hom_var.annotate(variant_category='hom_var')
merged_comphets = merged_comphets.annotate(variant_category='comphet')

merged_comphets_xlr_hom_var = merged_comphets.drop('proband_GT','proband_GT_set','proband_PBT_GT_set').union(xlr_phased).union(phased_hom_var)
# Annotate PAR status
merged_comphets_xlr_hom_var = merged_comphets_xlr_hom_var.annotate(in_non_par=~(merged_comphets_xlr_hom_var.locus.in_autosome_or_par()))
merged_comphets_xlr_hom_var = get_transmission(merged_comphets_xlr_hom_var)

merged_comphets_xlr_hom_var.flatten().export(f"{prefix}_{variant_types}_comp_hets_xlr_hom_var.tsv.gz")