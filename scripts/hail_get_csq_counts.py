import hail as hl

import gzip
import os
import glob
import io
import ast
import warnings
import seaborn as sns
import pandas as pd
import numpy as np
import sys

vcf_uri = sys.argv[1]
mpc_ht_uri = sys.argv[2]
ped_uri = sys.argv[3]
ad_threshold = int(sys.argv[4])
gq_threshold = int(sys.argv[5])
ab_min = float(sys.argv[6])
ab_max = float(sys.argv[7])
cores = sys.argv[8]  # string
mem = int(np.floor(float(sys.argv[9])))

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

mt = hl.import_vcf(vcf_uri, reference_genome='GRCh38', force_bgz=True, call_fields=[], array_elements_required=False)

mt = mt.annotate_entries(AB = mt.AD[1] / hl.sum(mt.AD))

mt = mt.filter_entries((mt.AD[1] >= ad_threshold) 
                       & (mt.GQ >= gq_threshold)
                       & (mt.AB >= ab_min) & (mt.AB <= ab_max))

ht = hl.read_table(mpc_ht_uri).key_by('locus', 'alleles')

mt = mt.annotate_rows(MPC=ht[mt.row_key].mpc)

# filter samples by pedigree
ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
ped.columns = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']

mt = mt.filter_cols(hl.array(ped['sample_id'].dropna().tolist()).contains(mt.s))

# annotate phenotype
ped_ht = hl.import_table(ped_uri).key_by('sample_id')
mt = mt.annotate_cols(phenotype=ped_ht[mt.s].phenotype)

# try this from https://github.com/Nealelab/recessive/blob/401af812dc4f1fc51cf2e8912aa598e2cef44e3c/Hail_%26_Export_Pipeline_Genotyped_dataset.ipynb
from typing import *

# Note that this is the current as of v81 with some included for backwards compatibility (VEP <= 75)
CSQ_CODING_HIGH_IMPACT = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost"]

CSQ_CODING_MEDIUM_IMPACT = [
    "start_lost",  # new in v81
    "initiator_codon_variant",  # deprecated
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",  # new in v79
    "splice_region_variant"
]

CSQ_CODING_LOW_IMPACT = [
    "incomplete_terminal_codon_variant",
    "start_retained_variant",  # new in v92
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant"]

CSQ_NON_CODING = [
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "non_coding_exon_variant",  # deprecated
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "nc_transcript_variant",  # deprecated
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "feature_elongation",
    "regulatory_region_variant",
    "feature_truncation",
    "intergenic_variant"
]

CSQ_ORDER = CSQ_CODING_HIGH_IMPACT + CSQ_CODING_MEDIUM_IMPACT + CSQ_CODING_LOW_IMPACT + CSQ_NON_CODING + ['.', 'NA']

def filter_vep_to_canonical_transcripts(mt: Union[hl.MatrixTable, hl.Table],
                                        vep_root: str = 'vep') -> Union[hl.MatrixTable, hl.Table]:
    canonical = mt[vep_root].transcript_consequences.filter(lambda csq: csq.CANONICAL == "YES" | csq.MANE_SELECT != '')
    vep_data = mt[vep_root].annotate(transcript_consequences=canonical)
    return mt.annotate_rows(**{vep_root: vep_data}) if isinstance(mt, hl.MatrixTable) else mt.annotate(**{vep_root: vep_data})

def add_most_severe_consequence_to_consequence(tc: hl.expr.StructExpression) -> hl.expr.StructExpression:

    """
    Add most_severe_consequence annotation to transcript consequences
    This is for a given transcript, as there are often multiple annotations for a single transcript:
    e.g. splice_region_variant&intron_variant -> splice_region_variant
    """

    csqs = hl.literal(CSQ_ORDER)

    return tc.annotate(
        most_severe_consequence=csqs.find(lambda c: tc.Consequence.contains(c))
    )

def process_consequences(mt: Union[hl.MatrixTable, hl.Table], vep_root: str = 'vep',
                         penalize_flags: bool = True) -> Union[hl.MatrixTable, hl.Table]:
    """
    Adds most_severe_consequence (worst consequence for a transcript) into [vep_root].transcript_consequences,
    and worst_csq_by_gene, any_LoF into [vep_root]
    :param MatrixTable mt: Input MT
    :param str vep_root: Root for vep annotation (probably vep)
    :param bool penalize_flags: Whether to penalize LoFTEE flagged variants, or treat them as equal to HC
    :return: MT with better formatted consequences
    :rtype: MatrixTable
    """
    csqs = hl.literal(CSQ_ORDER)
    csq_dict = hl.literal(dict(zip(CSQ_ORDER, range(len(CSQ_ORDER)))))

    def find_worst_transcript_consequence(tcl: hl.expr.ArrayExpression) -> hl.expr.StructExpression:
        """
        Gets worst transcript_consequence from an array of em
        """
        flag_score = 500
        no_flag_score = flag_score * (1 + penalize_flags)
        non_coding_score = 600
        non_canonical_score = 500
        def csq_score(tc):
            return csq_dict[csqs.find(lambda x: x == tc.most_severe_consequence)]
        
        # EDITED: PENALIZE NON-CODING AND NON-CANONICAL
        tcl = tcl.map(lambda tc: tc.annotate(
            csq_score=hl.case(missing_false=True)
            .when((tc.BIOTYPE != 'protein_coding'), csq_score(tc) + non_coding_score)
            .when((tc.CANONICAL != 'YES'), csq_score(tc) + non_canonical_score)
            .when((tc.LoF == 'HC') & (tc.LoF_flags == ''), csq_score(tc) - no_flag_score)
            .when((tc.LoF == 'HC') & (tc.LoF_flags != ''), csq_score(tc) - flag_score)
            .when(tc.LoF == 'LC', csq_score(tc) - 10) 
            .when(tc.PolyPhen.contains('probably_damaging'), csq_score(tc) - 0.5)  # EDITED
            .when(tc.PolyPhen.contains('possibly_damaging'), csq_score(tc) - 0.25)
            .when(tc.PolyPhen.contains('benign'), csq_score(tc) - 0.1)
            .default(csq_score(tc))
        ))
        return hl.or_missing(hl.len(tcl) > 0, hl.sorted(tcl, lambda x: x.csq_score)[0])

    transcript_csqs = mt[vep_root].transcript_consequences.map(add_most_severe_consequence_to_consequence)

    gene_dict = transcript_csqs.group_by(lambda tc: tc.SYMBOL)
    worst_csq_gene = gene_dict.map_values(find_worst_transcript_consequence).values()
    sorted_scores = hl.sorted(worst_csq_gene, key=lambda tc: tc.csq_score)
    lowest_score = hl.or_missing(hl.len(sorted_scores) > 0, sorted_scores[0].csq_score)
    gene_with_worst_csq = sorted_scores.filter(lambda tc: tc.csq_score == lowest_score).map(lambda tc: tc.SYMBOL)
    ensg_with_worst_csq = sorted_scores.filter(lambda tc: tc.csq_score == lowest_score).map(lambda tc: tc.Gene)

    vep_data = mt[vep_root].annotate(transcript_consequences=transcript_csqs,
                                     worst_consequence_term=csqs.find(lambda c: transcript_csqs.map(lambda csq: csq.most_severe_consequence).contains(c)),
                                     worst_csq_by_gene=sorted_scores,  # EDITED
                                     worst_csq=sorted_scores[0],
                                     any_LoF=hl.any(lambda x: x.LoF == 'HC', worst_csq_gene),
                                     gene_with_most_severe_csq=gene_with_worst_csq,
                                     ensg_with_most_severe_csq=ensg_with_worst_csq)
    
    return mt.annotate_rows(**{vep_root: vep_data}) if isinstance(mt, hl.MatrixTable) else mt.annotate(**{vep_root: vep_data})

header = hl.get_vcf_metadata(vcf_uri)
csq_columns = header['info']['CSQ']['Description'].split('Format: ')[1].split('|')

mt = mt.annotate_rows(vep=mt.info)
transcript_consequences = mt.vep.CSQ.map(lambda x: x.split('\|'))

transcript_consequences_strs = transcript_consequences.map(lambda x: hl.if_else(hl.len(x)>1, hl.struct(**
                                                       {col: x[i] if col!='Consequence' else x[i].split('&')  
                                                        for i, col in enumerate(csq_columns)}), 
                                                        hl.struct(**{col: hl.missing('str') if col!='Consequence' else hl.array([hl.missing('str')])  
                                                        for i, col in enumerate(csq_columns)})))

mt = mt.annotate_rows(vep=mt.vep.annotate(transcript_consequences=transcript_consequences_strs))
mt = mt.annotate_rows(vep=mt.vep.select('transcript_consequences'))

mt_csq = process_consequences(mt)

# misA = MPC >=1 and MPC <2
# MisB = MPC >=2

PTVs = ['frameshift_variant', 'stop_gained', 'splice_donor_variant', 'splice_acceptor_variant', 'transcript_ablation']
missense = ['missense_variant']
synonymous = ['synonymous_variant', 'stop_retained_variant']

new_csq_cats = {x: 'PTV' for x in PTVs} | {x: 'missense' for x in missense} | {x: 'synonymous' for x in synonymous}
mt_csq_filt = mt_csq.filter_rows(hl.array(PTVs + missense + synonymous).contains(mt_csq.vep.worst_consequence_term))

mt_csq_filt = mt_csq_filt.annotate_rows(overall_csq_term = hl.dict(new_csq_cats)[mt_csq_filt.vep.worst_consequence_term])
mt_csq_filt = mt_csq_filt.annotate_rows(overall_csq_term = hl.if_else(mt_csq_filt.overall_csq_term=='missense',
                                                                     hl.case()
                                                                      .when((mt_csq_filt.MPC>=1)&(mt_csq_filt.MPC<2), 'misA')
                                                                      .when(mt_csq_filt.MPC>=2, 'misB').or_missing(),
                                                                     mt_csq_filt.overall_csq_term))

mt_csq_filt = mt_csq_filt.filter_rows(hl.is_defined(mt_csq_filt.overall_csq_term))

mt_csq_filt_controls = mt_csq_filt.filter_cols(mt_csq_filt.phenotype=='1')
mt_csq_filt_controls = hl.variant_qc(mt_csq_filt_controls)
mt_csq_filt_controls = mt_csq_filt_controls.filter_rows(mt_csq_filt_controls.variant_qc.AC[1]>0)
mt_csq_controls_rows = mt_csq_filt_controls.rows()

count_mt_controls = mt_csq_controls_rows.group_by(mt_csq_controls_rows.vep.worst_csq.SYMBOL)\
.aggregate(csq_count=hl.agg.counter(mt_csq_controls_rows.overall_csq_term))

count_df_controls = count_mt_controls.to_pandas()

mt_csq_filt_cases = mt_csq_filt.filter_cols(mt_csq_filt.phenotype=='2')
mt_csq_filt_cases = hl.variant_qc(mt_csq_filt_cases)
mt_csq_filt_cases = mt_csq_filt_cases.filter_rows(mt_csq_filt_cases.variant_qc.AC[1]>0)
mt_csq_cases_rows = mt_csq_filt_cases.rows()

count_mt_cases = mt_csq_cases_rows.group_by(mt_csq_cases_rows.vep.worst_csq.SYMBOL)\
.aggregate(csq_count=hl.agg.counter(mt_csq_cases_rows.overall_csq_term))

count_df_cases = count_mt_cases.to_pandas()

count_df_controls_tmp = pd.json_normalize(count_df_controls.csq_count)
count_df_controls_tmp.columns = count_df_controls_tmp.columns + '_controls'
controls_tmp = pd.concat([count_df_controls.SYMBOL, count_df_controls_tmp], axis=1)

count_df_cases_tmp = pd.json_normalize(count_df_cases.csq_count)
count_df_cases_tmp.columns = count_df_cases_tmp.columns + '_cases'
cases_tmp = pd.concat([count_df_cases.SYMBOL, count_df_cases_tmp], axis=1)

final_counts = controls_tmp.merge(cases_tmp)

final_counts.to_csv(os.path.basename(vcf_uri).split('.vcf')[0] + '_mis_syn_ptv_counts.tsv', sep='\t', index=False)