version 1.0

import "wes-denovo-helpers.wdl" as helpers
import "scatterHailMTs.wdl" as scatterHailMTs

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow exportVDStoVCF {
    input {
        File sample_file
        String input_vds
        String hail_docker
        String output_vcf_basename
        Int n_shards
    }
    
    call helpers.getHailMTSize as getInputVDSSize {
        input:
            mt_uri=input_vds,
            hail_docker=hail_docker
    }

    call scatterHailMTs.getRepartitions as getRepartitions {
        input:
        n_shards=n_shards,
        mt_uri=input_vds,
        hail_docker=hail_docker
    }

    scatter (interval in getRepartitions.partition_intervals) {
        call exportVDS {
            input:
                sample_file=sample_file,
                input_vds=input_vds,
                output_vcf_basename=output_vcf_basename,
                shard_n=interval[0],
                interval_start=interval[1],
                interval_end=interval[2],
                hail_docker=hail_docker,
                input_size=getInputVDSSize.mt_size / n_shards
        }
    }


    output {
        Array[File] vcf_shards = exportVDS.vcf_shard
        Array[File] vcf_shards_index = exportVDS.vcf_shard_idx
    }
}

task exportVDS {
    input {
        File sample_file
        String input_vds
        String output_vcf_basename
        Int shard_n
        Int interval_start
        Int interval_end
        String hail_docker
        Float input_size
        RuntimeAttr? runtime_attr_override
    }

    Float base_disk_gb = 10.0
    Float input_disk_scale = 2.0

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
    cat <<EOF > export_vds.py
    import pandas as pd
    import numpy as np
    import hail as hl
    import os
    import sys

    input_vds = sys.argv[1]
    output_vcf_basename = sys.argv[2]
    sample_file = sys.argv[3]
    shard_n = int(sys.argv[4])
    interval_start = int(sys.argv[5])
    interval_end = int(sys.argv[6])
    cores = sys.argv[7]  # string
    mem = int(np.floor(float(sys.argv[8])))

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

    # move gvcf_info from entries to rows before merging VDS
    variant_mt = vds.variant_data
    rows = variant_mt.entries().select('rsid','gvcf_info').key_by('locus', 'alleles')
    variant_mt = variant_mt.annotate_rows(info=rows[variant_mt.row_key].gvcf_info).drop('gvcf_info')
    vds.variant_data = variant_mt

    # merge VDS with specific range of partitions
    mt = to_dense_mt(vds, interval_start, interval_end)

    # subset samples
    mt = mt.filter_cols(hl.array(samples).contains(mt.s))
    # convert LGT to GT
    mt = mt.annotate_entries(GT=hl.vds.lgt_to_gt(mt.LGT, mt.LA))

    # remove all AC=0
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.AC[1] > 0, keep = True)
    mt = mt.annotate_rows(info=mt.info.annotate(AC=mt.variant_qc.AC[1:], 
                                       AF=mt.variant_qc.AF[1:],
                                       AN=mt.variant_qc.AN))
    mt = mt.drop('variant_qc')

    hl.export_vcf(mt, f"{output_vcf_basename}_shard_{shard_n}.vcf.bgz", tabix=True)
    EOF
    python3 export_vds.py ~{input_vds} ~{output_vcf_basename} ~{sample_file} \
        ~{shard_n} ~{interval_start} ~{interval_end} ~{cpu_cores} ~{memory}
    >>>

    output {
        File vcf_shard = "~{output_vcf_basename}_shard_~{shard_n}.vcf.bgz"
        File vcf_shard_idx = "~{output_vcf_basename}_shard_~{shard_n}.vcf.bgz.tbi"
    }
}
