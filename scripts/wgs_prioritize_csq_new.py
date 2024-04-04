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

hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

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


def process_consequence_cohort(csq_columns, vcf_metrics_uri, numeric, sample_column):
    ht = hl.import_table(vcf_metrics_uri)

    ht = ht.annotate(locus=hl.parse_variant(ht.ID, reference_genome='GRCh38').locus,
                            alleles=hl.parse_variant(ht.ID, reference_genome='GRCh38').alleles)

    mt = ht.to_matrix_table_row_major(columns=numeric, entry_field_name='metrics', col_field_name=sample_column)
    mt = mt.key_rows_by(mt.locus, mt.alleles)

    transcript_consequences = mt.CSQ.replace("\[",'').replace("\]",'').replace("\'",'').replace(' ','').split(',').map(lambda x: x.split('\|'))

    transcript_consequences_strs = transcript_consequences.map(lambda x: hl.if_else(hl.len(x)>1, hl.struct(**
                                                           {col: x[i] if col!='Consequence' else x[i].split('&')  
                                                            for i, col in enumerate(csq_columns)}), 
                                                            hl.struct(**{col: '.' if col!='Consequence' else hl.array(['.'])  
                                                            for i, col in enumerate(csq_columns)})))

    mt_filt=mt.annotate_rows(vep=hl.Struct(transcript_consequences = transcript_consequences_strs))

    # mt_filt = filter_vep_to_canonical_transcripts(mt)
    mt_csq = process_consequences(mt_filt)
    mt_csq_annotated = mt_csq.annotate_rows(worst_csq=mt_csq.vep.worst_csq)
    mt_csq_annotated = mt_csq_annotated.drop('vep')

    mt_csq_annotated = mt_csq_annotated.annotate_rows(
        isPTV = mt_csq_annotated.worst_csq.Consequence.contains('frameshift_variant') | mt_csq_annotated.worst_csq.Consequence.contains('stop_gained') |
                mt_csq_annotated.worst_csq.Consequence.contains('splice_donor_variant') | mt_csq_annotated.worst_csq.Consequence.contains('splice_acceptor_variant') |
                mt_csq_annotated.worst_csq.Consequence.contains('transcript_ablation'),
        isMIS = mt_csq_annotated.worst_csq.Consequence.contains('missense_variant'),
        isSYN = mt_csq_annotated.worst_csq.Consequence.contains('synonymous_variant') | mt_csq_annotated.worst_csq.Consequence.contains('stop_retained_variant'),
        isSRG = mt_csq_annotated.worst_csq.Consequence.contains('splice_region_variant'),
        isSSL = mt_csq_annotated.worst_csq.Consequence.contains('start_lost') | mt_csq_annotated.worst_csq.Consequence.contains('stop_lost'),
        isINF = mt_csq_annotated.worst_csq.Consequence.contains('inframe_insertion') | mt_csq_annotated.worst_csq.Consequence.contains('inframe_deletion') )

    # And variants have an interesting "most severe" consequence (SII = "severe is interesting")
    mt_csq_annotated = mt_csq_annotated.annotate_rows(
        SII = (mt_csq_annotated.worst_csq.most_severe_consequence == 'frameshift_variant') | (mt_csq_annotated.worst_csq.most_severe_consequence == 'stop_gained') |
            (mt_csq_annotated.worst_csq.most_severe_consequence == 'splice_donor_variant') | (mt_csq_annotated.worst_csq.most_severe_consequence ==  'splice_acceptor_variant') |
            (mt_csq_annotated.worst_csq.most_severe_consequence == 'transcript_ablation') | (mt_csq_annotated.worst_csq.most_severe_consequence == 'missense_variant') |
            (mt_csq_annotated.worst_csq.most_severe_consequence == 'synonymous_variant') | (mt_csq_annotated.worst_csq.most_severe_consequence == 'stop_retained_variant') |
            (mt_csq_annotated.worst_csq.most_severe_consequence == 'splice_region_variant') | (mt_csq_annotated.worst_csq.most_severe_consequence == 'start_lost') | 
            (mt_csq_annotated.worst_csq.most_severe_consequence == 'stop_lost') | (mt_csq_annotated.worst_csq.most_severe_consequence == 'inframe_insertion') | 
            (mt_csq_annotated.worst_csq.most_severe_consequence == 'inframe_deletion') )

    # Additionally annotate SNV vs indel
    mt_csq_annotated = mt_csq_annotated.annotate_rows(isSNV = hl.is_snp(mt_csq_annotated.alleles[0], mt_csq_annotated.alleles[1]),
                        isIndel = hl.is_indel(mt_csq_annotated.alleles[0], mt_csq_annotated.alleles[1]))

    mt_csq_rows = mt_csq_annotated.rows().flatten()
    mt_csq_df = mt_csq_rows.to_pandas()

    vcf_metrics = mt_csq_df.rename({col: col.split('worst_csq.')[1] for col in mt_csq_df.columns if 'worst_csq' in col}, axis=1)
    return vcf_metrics

csq_columns_less = ['Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 
                    'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 
                    'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'ALLELE_NUM', 
                    'DISTANCE', 'STRAND', 'FLAGS', 'VARIANT_CLASS', 'MINIMISED', 'SYMBOL_SOURCE', 
                    'HGNC_ID', 'CANONICAL', 'TSL', 'APPRIS', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 
                    'UNIPARC', 'GENE_PHENO', 'SIFT', 'PolyPhen', 'DOMAINS', 'miRNA', 'HGVS_OFFSET', 
                    'AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'AA_AF', 'EA_AF', 'gnomAD_AF', 
                    'gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 
                    'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF', 'MAX_AF', 'MAX_AF_POPS', 'CLIN_SIG', 
                    'SOMATIC', 'PHENO', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 
                    'LOEUF', 'LoF', 'LoF_filter', 'LoF_flags', 'LoF_info']

csq_columns_more = ["Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature",
                   "BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position",
                   "Protein_position","Amino_acids","Codons","Existing_variation","ALLELE_NUM",
                   "DISTANCE","STRAND","FLAGS","VARIANT_CLASS","MINIMISED","SYMBOL_SOURCE",
                   "HGNC_ID","CANONICAL","MANE_SELECT","MANE_PLUS_CLINICAL","TSL","APPRIS",
                   "CCDS","ENSP","SWISSPROT","TREMBL","UNIPARC","UNIPROT_ISOFORM","GENE_PHENO",
                   "SIFT","PolyPhen","DOMAINS","miRNA","HGVS_OFFSET","AF","AFR_AF","AMR_AF",
                   "EAS_AF","EUR_AF","SAS_AF","gnomADe_AF","gnomADe_AFR_AF","gnomADe_AMR_AF",
                   "gnomADe_ASJ_AF","gnomADe_EAS_AF","gnomADe_FIN_AF","gnomADe_NFE_AF","gnomADe_OTH_AF",
                   "gnomADe_SAS_AF","gnomADg_AF","gnomADg_AFR_AF","gnomADg_AMI_AF","gnomADg_AMR_AF",
                   "gnomADg_ASJ_AF","gnomADg_EAS_AF","gnomADg_FIN_AF","gnomADg_MID_AF","gnomADg_NFE_AF",
                   "gnomADg_OTH_AF","gnomADg_SAS_AF","MAX_AF","MAX_AF_POPS","CLIN_SIG","SOMATIC",
                   "PHENO","PUBMED","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE",
                   "TRANSCRIPTION_FACTORS","LoF","LoF_filter","LoF_flags","LoF_info"]

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

final_output = pd.read_csv(vcf_metrics_uri, sep='\t')
final_output['CSQ'] = final_output.CSQ.replace({'.':np.nan}).str.split(',')
n_csq_fields = len(final_output[~final_output.CSQ.isna()].CSQ.iloc[0][0].split('|'))

if n_csq_fields==len(csq_columns_more):
    csq_columns = csq_columns_more
elif n_csq_fields==len(csq_columns_less):
    csq_columns = csq_columns_less
else:
    warnings.simplefilter("error")
    warnings.warn("CSQ fields are messed up!")

numeric = []
for col in final_output.columns:
    if col==sample_column:
        continue
    try:
        final_output[col].astype(float)
        numeric.append(col)
    except:
        continue

df = process_consequence_cohort(csq_columns, vcf_metrics_uri, numeric, sample_column)
df['isCoding'] = df.Consequence.astype(str).replace({'None': '[]'}).apply(ast.literal_eval).apply(lambda csq: np.intersect1d(csq, coding_variants).size!=0)

df = pd.concat([final_output, df[np.setdiff1d(df.columns, final_output.columns)]], axis=1)
df.to_csv(f"{os.path.basename(vcf_metrics_uri).split('.tsv')[0]}_prioritized_csq.tsv", sep='\t',index=False)