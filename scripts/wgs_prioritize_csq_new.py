import hail as hl
import pandas as pd
import numpy as np
import sys
import os
import warnings
import ast

vcf_metrics_uri = sys.argv[1]
cores = sys.argv[2]
mem = int(np.floor(float(sys.argv[3])))
sample_column = sys.argv[4]
vep_vcf_uri = sys.argv[5]
genome_build = sys.argv[6]

hl.init(min_block_size=128, 
        local=f"local[*]", 
        spark_conf={
                    "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.executor.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.speculation": 'true'}, 
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

# Prioritized CSQ

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
    canonical = mt[vep_root].transcript_consequences.filter(lambda csq: csq.CANONICAL == "YES")
    vep_data = mt[vep_root].annotate(transcript_consequences=hl.if_else(canonical.size()>0, canonical, mt[vep_root].transcript_consequences))
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
            # .when((tc.LoF == 'HC') & (tc.LoF_flags == ''), csq_score(tc) - no_flag_score)
            # .when((tc.LoF == 'HC') & (tc.LoF_flags != ''), csq_score(tc) - flag_score)
            # .when(tc.LoF == 'LC', csq_score(tc) - 10) 
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
                                    #  any_LoF=hl.any(lambda x: x.LoF == 'HC', worst_csq_gene),
                                     gene_with_most_severe_csq=gene_with_worst_csq,
                                     ensg_with_most_severe_csq=ensg_with_worst_csq)
    
    return mt.annotate_rows(**{vep_root: vep_data}) if isinstance(mt, hl.MatrixTable) else mt.annotate(**{vep_root: vep_data})

coding_variants = ['coding_sequence_variant', 'frameshift_variant', 
        'incomplete_terminal_codon_variant', 'inframe_deletion', 'inframe_insertion',
        'missense_variant', 'protein_altering_variant', 'splice_acceptor_variant',
        'splice_donor_variant', 'start_lost', 'stop_gained', 'stop_lost',
        'stop_retained_variant', 'synonymous_variant']

noncoding_variants = ['3_prime_UTR_variant', '5_prime_UTR_variant',
        'downstream_gene_variant', 'intergenic_variant', 'intron_variant', 
        'mature_miRNA_variant', 'non_coding_transcript_exon_variant',
        'splice_region_variant', 'upstream_gene_variant']

PTVs = ['frameshift_variant', 'stop_gained', 'splice_donor_variant', 'splice_acceptor_variant', 'transcript_ablation']
missense = ['missense_variant']
synonymous = ['synonymous_variant', 'stop_retained_variant']

final_output = pd.read_csv(vcf_metrics_uri, sep='\t', compression="gzip" if vcf_metrics_uri.split('.')[-1] in ['gz', 'bgz'] else None)
# force_bgz broken in hl.import_table, so need to uncompress
final_output.to_csv(f"{os.path.basename(vcf_metrics_uri).split('.tsv')[0]}_uncompressed.tsv", sep='\t', index=False)
final_output['CSQ'] = final_output.CSQ.replace({'.':np.nan}).str.split(',')

header = hl.get_vcf_metadata(vep_vcf_uri)
csq_columns = header['info']['CSQ']['Description'].split('Format: ')[1].split('|')

ht = hl.import_table(f"{os.path.basename(vcf_metrics_uri).split('.tsv')[0]}_uncompressed.tsv")

if 'ID' not in list(ht.row):
    ht = ht.annotate(alleles=ht.alleles.replace("\[",'').replace("\]",'').replace("\'",'').replace(' ','').replace('"','').split(','))
    ht = ht.annotate(ID=ht.locus+':'+ht.alleles[0]+':'+ht.alleles[1])

ht = ht.annotate(locus=hl.parse_variant(ht.ID, reference_genome=genome_build).locus,
                        alleles=hl.parse_variant(ht.ID, reference_genome=genome_build).alleles)

if 'VarKey' not in list(ht.row):
    ht = ht.annotate(VarKey=ht.ID+':'+ht[sample_column])

ht = ht.select('locus','alleles','VarKey','CSQ')
transcript_consequences = ht.CSQ.replace("\[",'').replace("\]",'').replace("\'",'').replace(' ','').replace('"','').split(',').map(lambda x: x.split('\|'))

transcript_consequences_strs = transcript_consequences.map(lambda x: hl.if_else(hl.len(x)>1, hl.struct(**
                                                       {col: x[i] if col!='Consequence' else x[i].split('&')  
                                                        for i, col in enumerate(csq_columns)}), 
                                                        hl.struct(**{col: '.' if col!='Consequence' else hl.array(['.'])  
                                                        for i, col in enumerate(csq_columns)})))

ht=ht.annotate(vep=hl.Struct(transcript_consequences = transcript_consequences_strs))

# filter canonical
ht = filter_vep_to_canonical_transcripts(ht)

ht_csq = process_consequences(ht)
ht_csq = ht_csq.annotate(worst_csq=ht_csq.vep.worst_csq)
ht_csq = ht_csq.drop('vep')

ht_csq = ht_csq.annotate(
    isPTV = ht_csq.worst_csq.Consequence.contains('frameshift_variant') | ht_csq.worst_csq.Consequence.contains('stop_gained') |
            ht_csq.worst_csq.Consequence.contains('splice_donor_variant') | ht_csq.worst_csq.Consequence.contains('splice_acceptor_variant') |
            ht_csq.worst_csq.Consequence.contains('transcript_ablation'),
    isMIS = ht_csq.worst_csq.Consequence.contains('missense_variant'),
    isSYN = ht_csq.worst_csq.Consequence.contains('synonymous_variant') | ht_csq.worst_csq.Consequence.contains('stop_retained_variant'),
    isSRG = ht_csq.worst_csq.Consequence.contains('splice_region_variant'),
    isSSL = ht_csq.worst_csq.Consequence.contains('start_lost') | ht_csq.worst_csq.Consequence.contains('stop_lost'),
    isINF = ht_csq.worst_csq.Consequence.contains('inframe_insertion') | ht_csq.worst_csq.Consequence.contains('inframe_deletion') )

# And variants have an interesting "most severe" consequence (SII = "severe is interesting")
ht_csq = ht_csq.annotate(
    SII = (ht_csq.worst_csq.most_severe_consequence == 'frameshift_variant') | (ht_csq.worst_csq.most_severe_consequence == 'stop_gained') |
        (ht_csq.worst_csq.most_severe_consequence == 'splice_donor_variant') | (ht_csq.worst_csq.most_severe_consequence ==  'splice_acceptor_variant') |
        (ht_csq.worst_csq.most_severe_consequence == 'transcript_ablation') | (ht_csq.worst_csq.most_severe_consequence == 'missense_variant') |
        (ht_csq.worst_csq.most_severe_consequence == 'synonymous_variant') | (ht_csq.worst_csq.most_severe_consequence == 'stop_retained_variant') |
        (ht_csq.worst_csq.most_severe_consequence == 'splice_region_variant') | (ht_csq.worst_csq.most_severe_consequence == 'start_lost') | 
        (ht_csq.worst_csq.most_severe_consequence == 'stop_lost') | (ht_csq.worst_csq.most_severe_consequence == 'inframe_insertion') | 
        (ht_csq.worst_csq.most_severe_consequence == 'inframe_deletion') )

# Additionally annotate SNV vs indel
ht_csq = ht_csq.annotate(isSNV = hl.is_snp(ht_csq.alleles[0], ht_csq.alleles[1]),
                    isIndel = hl.is_indel(ht_csq.alleles[0], ht_csq.alleles[1]))
ht_csq.flatten().export(f"{os.path.basename(vcf_metrics_uri).split('.tsv')[0]}_temp.tsv.gz")

ht_csq_df = pd.read_csv(f"{os.path.basename(vcf_metrics_uri).split('.tsv')[0]}_temp.tsv.gz", sep='\t')
df = ht_csq_df.rename({col: col.split('worst_csq.')[1] for col in ht_csq_df.columns if 'worst_csq' in col}, axis=1)
df.index = df['VarKey']

df['isCoding'] = df.Consequence.astype(str).replace({'None': '[]'}).apply(ast.literal_eval).apply(lambda csq: np.intersect1d(csq, coding_variants).size!=0)

if 'ID' not in final_output.columns:
    final_output['ID'] = final_output.locus.astype(str) + ':' + final_output.alleles.apply(ast.literal_eval).apply(':'.join)
if 'VarKey' not in final_output.columns:
    final_output['VarKey'] = final_output[['ID', sample_column]].astype(str).agg(':'.join, axis=1)

final_output.index = final_output['VarKey']

new_cols = list(np.setdiff1d(df.columns, ['locus','alleles','VarKey','CSQ']))
final_output[new_cols] = df[new_cols]  # allows overwriting columns

final_output.to_csv(f"{os.path.basename(vcf_metrics_uri).split('.tsv')[0]}_prioritized_csq.tsv.gz", sep='\t',index=False)