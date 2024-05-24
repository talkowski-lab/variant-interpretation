import pandas as pd
import numpy as np
import hail as hl
import os
import sys

input_vds = sys.argv[1]
output_vcf_basename = sys.argv[2]
sample_file = sys.argv[3]
info_ht_uri = sys.argv[4]
vep_ht_uri = sys.argv[5]
qc_ht_uri = sys.argv[6]
shard_n = int(sys.argv[7])
interval_start = int(sys.argv[8])
interval_end = int(sys.argv[9])
cores = sys.argv[10]  # string
mem = int(np.floor(float(sys.argv[11])))

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

samples = pd.read_csv(sample_file, header=None)[0].tolist()

vds = hl.vds.read_vds(input_vds)
from hail import ir
from hail.expr import expr_any, expr_array, expr_bool, expr_interval, expr_locus, expr_str
from hail.matrixtable import MatrixTable
from hail.table import Table
from hail.typecheck import dictof, enumeration, func_spec, nullable, oneof, sequenceof, typecheck
from hail.utils.java import Env, info, warning
from hail.utils.misc import new_temp_file, wrap_to_list
from hail.vds.variant_dataset import VariantDataset

# FUNC START
@typecheck(vds=VariantDataset, interval_start=int, interval_end=int)
def to_dense_mt(vds: 'VariantDataset', interval_start: int, interval_end: int) -> 'MatrixTable':
    """Creates a single, dense :class:`.MatrixTable` from the split
    :class:`.VariantDataset` representation.

    Parameters
    ----------
    vds : :class:`.VariantDataset`
        Dataset in VariantDataset representation.

    Returns
    -------
    :class:`.MatrixTable`
        Dataset in dense MatrixTable representation.
    """
    ref = vds.reference_data._filter_partitions((interval_start, interval_end))  # EDITED
    # FIXME(chrisvittal) consider changing END semantics on VDS to make this better
    # see https://github.com/hail-is/hail/issues/13183 for why this is here and more discussion
    # we assume that END <= contig.length
    ref = ref.annotate_rows(_locus_global_pos=ref.locus.global_position(), _locus_pos=ref.locus.position)
    ref = ref.transmute_entries(_END_GLOBAL=ref._locus_global_pos + (ref.END - ref._locus_pos))

    to_drop = 'alleles', 'rsid', 'ref_allele', '_locus_global_pos', '_locus_pos'
    ref = ref.drop(*(x for x in to_drop if x in ref.row))
    var = vds.variant_data._filter_partitions((interval_start, interval_end))  # EDITED
    refl = ref.localize_entries('_ref_entries')
    varl = var.localize_entries('_var_entries', '_var_cols')
    varl = varl.annotate(_variant_defined=True)
    joined = varl.key_by('locus').join(refl, how='outer')
    dr = joined.annotate(
        dense_ref=hl.or_missing(
            joined._variant_defined, hl.scan._densify(hl.len(joined._var_cols), joined._ref_entries)
        )
    )
    dr = dr.filter(dr._variant_defined)

    def coalesce_join(ref, var):
        call_field = 'GT' if 'GT' in var else 'LGT'
        assert call_field in var, var.dtype

        shared_fields = [call_field, *list(f for f in ref.dtype if f in var.dtype)]
        shared_field_set = set(shared_fields)
        var_fields = [f for f in var.dtype if f not in shared_field_set]

        return hl.if_else(
            hl.is_defined(var),
            var.select(*shared_fields, *var_fields),
            ref.annotate(**{call_field: hl.call(0, 0)}).select(
                *shared_fields, **{f: hl.missing(var[f].dtype) for f in var_fields}
            ),
        )

    dr = dr.annotate(
        _dense=hl.rbind(
            dr._ref_entries,
            lambda refs_at_this_row: hl.enumerate(hl.zip(dr._var_entries, dr.dense_ref)).map(
                lambda tup: coalesce_join(
                    hl.coalesce(
                        refs_at_this_row[tup[0]],
                        hl.or_missing(tup[1][1]._END_GLOBAL >= dr.locus.global_position(), tup[1][1]),
                    ),
                    tup[1][0],
                )
            ),
        ),
    )

    dr = dr._key_by_assert_sorted('locus', 'alleles')
    fields_to_drop = ['_var_entries', '_ref_entries', 'dense_ref', '_variant_defined']

    if hl.vds.VariantDataset.ref_block_max_length_field in dr.globals:
        fields_to_drop.append(hl.vds.VariantDataset.ref_block_max_length_field)

    if 'ref_allele' in dr.row:
        fields_to_drop.append('ref_allele')
    dr = dr.drop(*fields_to_drop)
    return dr._unlocalize_entries('_dense', '_var_cols', list(var.col_key))
# FUNC END

# merge VDS with specific range of partitions
mt = to_dense_mt(vds, interval_start, interval_end)

# subset samples
mt = mt.filter_cols(hl.array(samples).contains(mt.s))
# convert LGT to GT
mt = mt.annotate_entries(GT=hl.vds.lgt_to_gt(mt.LGT, mt.LA))
# split multi-allelic
mt = hl.split_multi(mt)

# get row/INFO fields from INFO HT
info_ht = hl.read_table(info_ht_uri)
mt = mt.annotate_rows(info=info_ht[mt.row_key].info)
mt = mt.drop('gvcf_info')

# get QUAL/FILTER info from QC HT
qc_ht = hl.read_table(qc_ht_uri)
mt = mt.annotate_rows(qual=qc_ht[mt.row_key].qual, filters=qc_ht[mt.row_key].filters)

# get VEP info
vep_ht = hl.read_table(vep_ht_uri)
mt = mt.annotate_rows(info=mt.info.annotate(vep=vep_ht[mt.row_key].vep))

# remove all AC=0
mt = hl.variant_qc(mt)
mt = mt.filter_rows(mt.variant_qc.AC[1] > 0, keep = True)
mt = mt.annotate_rows(info=mt.info.annotate(AC=mt.variant_qc.AC[1:], 
                                    AF=mt.variant_qc.AF[1:],
                                    AN=mt.variant_qc.AN))
mt = mt.drop('variant_qc')

hl.export_vcf(mt, f"{output_vcf_basename}_shard_{shard_n}.vcf.bgz", tabix=True)
