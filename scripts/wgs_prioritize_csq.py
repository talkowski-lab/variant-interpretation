import hail as hl
import pandas as pd
import numpy as np
import sys
import os
import warnings

vcf_metrics_uri = sys.argv[1]
cores = sys.argv[2]
mem = int(np.floor(float(sys.argv[3])))

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

CSQ_ORDER = CSQ_CODING_HIGH_IMPACT + CSQ_CODING_MEDIUM_IMPACT + CSQ_CODING_LOW_IMPACT + CSQ_NON_CODING

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

        def csq_score(tc):
            return csq_dict[csqs.find(lambda x: x == tc.most_severe_consequence)]
        tcl = tcl.map(lambda tc: tc.annotate(
            csq_score=hl.case(missing_false=True)
            .when((tc.LoF == 'HC') & (tc.LoF_flags == ''), csq_score(tc) - no_flag_score)
            .when((tc.LoF == 'HC') & (tc.LoF_flags != ''), csq_score(tc) - flag_score)
            .when(tc.LoF == 'LC', csq_score(tc) - 10)
            .when(tc.PolyPhen == 'probably_damaging', csq_score(tc) - 0.5)
            .when(tc.PolyPhen == 'possibly_damaging', csq_score(tc) - 0.25)
            .when(tc.PolyPhen == 'benign', csq_score(tc) - 0.1)
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
                                     worst_csq_by_gene=worst_csq_gene,
                                     any_LoF=hl.any(lambda x: x.LoF == 'HC', worst_csq_gene),
                                     gene_with_most_severe_csq=gene_with_worst_csq,
                                     ensg_with_most_severe_csq=ensg_with_worst_csq)

    return mt.annotate_rows(**{vep_root: vep_data}) if isinstance(mt, hl.MatrixTable) else mt.annotate(**{vep_root: vep_data})

#define VEP annotations
def add_vep_annotations(mt):
    
    # First, annotate based on parsing VEP output
    mt = mt.annotate_rows(
        SYMBOL = hl.coalesce(
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") & (hl.is_defined(x.Amino_acids)) ).SYMBOL,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") ).SYMBOL,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") ).SYMBOL,
            mt.vep.transcript_consequences.find(lambda x: x.Consequence.contains(mt.vep.worst_consequence_term)).SYMBOL),
        Gene = hl.coalesce(
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") & (hl.is_defined(x.Amino_acids)) ).Gene,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") ).Gene,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") ).Gene,
            mt.vep.transcript_consequences.find(lambda x: x.Consequence.contains(mt.vep.worst_consequence_term)).Gene),
        transcript_id = hl.coalesce(
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") & (hl.is_defined(x.Amino_acids)) ).Feature,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") ).Feature,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") ).Feature,
            mt.vep.transcript_consequences.find(lambda x: x.Consequence.contains(mt.vep.worst_consequence_term)).Feature),
        HGNC_ID = hl.coalesce(
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") & (hl.is_defined(x.Amino_acids)) ).HGNC_ID,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") ).HGNC_ID,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") ).HGNC_ID,
            mt.vep.transcript_consequences.find(lambda x: x.Consequence.contains(mt.vep.worst_consequence_term)).HGNC_ID),
        HGVSc = hl.coalesce(
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") & (hl.is_defined(x.Amino_acids)) ).HGVSc,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") ).HGVSc,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") ).HGVSc,
            mt.vep.transcript_consequences.find(lambda x: x.Consequence.contains(mt.vep.worst_consequence_term)).HGVSc),
        HGVSp = hl.coalesce(
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") & (hl.is_defined(x.Amino_acids)) ).HGVSp,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") ).HGVSp,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") ).HGVSp,
            mt.vep.transcript_consequences.find(lambda x: x.Consequence.contains(mt.vep.worst_consequence_term)).HGVSp),
        Consequence = hl.coalesce(
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") & (hl.is_defined(x.Amino_acids)) ).Consequence,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") ).Consequence,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") ).Consequence,
            mt.vep.transcript_consequences.find(lambda x: x.Consequence.contains(mt.vep.worst_consequence_term)).Consequence),
        MAX_AF = hl.coalesce(
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") & (hl.is_defined(x.Amino_acids)) ).MAX_AF,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") ).MAX_AF,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") ).MAX_AF,
            mt.vep.transcript_consequences.find(lambda x: x.Consequence.contains(mt.vep.worst_consequence_term)).MAX_AF),
        LoF = hl.coalesce(
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") & (hl.is_defined(x.Amino_acids)) ).LoF,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") ).LoF,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") ).LoF,
            mt.vep.transcript_consequences.find(lambda x: x.Consequence.contains(mt.vep.worst_consequence_term)).LoF),
        LoF_flags = hl.coalesce(
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") & (hl.is_defined(x.Amino_acids)) ).LoF_flags,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") ).LoF_flags,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") ).LoF_flags,
            mt.vep.transcript_consequences.find(lambda x: x.Consequence.contains(mt.vep.worst_consequence_term)).LoF_flags),
        PolyPhen = hl.coalesce(
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") & (hl.is_defined(x.Amino_acids)) ).PolyPhen,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") ).PolyPhen,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") ).PolyPhen,
            mt.vep.transcript_consequences.find(lambda x: x.Consequence.contains(mt.vep.worst_consequence_term)).PolyPhen),
        SIFT = hl.coalesce(
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") & (hl.is_defined(x.Amino_acids)) ).SIFT,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") & (x.BIOTYPE == "protein_coding") ).SIFT,
            mt.vep.transcript_consequences.find(lambda x: (x.CANONICAL == "YES") ).SIFT,
            mt.vep.transcript_consequences.find(lambda x: x.Consequence.contains(mt.vep.worst_consequence_term)).SIFT) )

    # Now, annotate based on the "consequence" annotation
    mt = mt.annotate_rows(
        isPTV = mt.Consequence.contains('frameshift_variant') | mt.Consequence.contains('stop_gained') |
                mt.Consequence.contains('splice_donor_variant') | mt.Consequence.contains('splice_acceptor_variant') |
                mt.Consequence.contains('transcript_ablation'),
        isMIS = mt.Consequence.contains('missense_variant'),
        isSYN = mt.Consequence.contains('synonymous_variant') | mt.Consequence.contains('stop_retained_variant'),
        isSRG = mt.Consequence.contains('splice_region_variant'),
        isSSL = mt.Consequence.contains('start_lost') | mt.Consequence.contains('stop_lost'),
        isINF = mt.Consequence.contains('inframe_insertion') | mt.Consequence.contains('inframe_deletion') )
    
    # And variants have an interesting "most severe" consequence (SII = "severe is interesting")
    mt = mt.annotate_rows(
        SII = (mt.vep.worst_consequence_term == 'frameshift_variant') | (mt.vep.worst_consequence_term == 'stop_gained') |
              (mt.vep.worst_consequence_term == 'splice_donor_variant') | (mt.vep.worst_consequence_term ==  'splice_acceptor_variant') |
              (mt.vep.worst_consequence_term == 'transcript_ablation') | (mt.vep.worst_consequence_term == 'missense_variant') |
              (mt.vep.worst_consequence_term == 'synonymous_variant') | (mt.vep.worst_consequence_term == 'stop_retained_variant') |
              (mt.vep.worst_consequence_term == 'splice_region_variant') | (mt.vep.worst_consequence_term == 'start_lost') | 
              (mt.vep.worst_consequence_term == 'stop_lost') | (mt.vep.worst_consequence_term == 'inframe_insertion') | 
              (mt.vep.worst_consequence_term == 'inframe_deletion') )
    
    # Additionally annotate SNP vs indel
    mt = mt.annotate_rows(isSNP = hl.is_snp(mt.alleles[0], mt.alleles[1]),
                          isIndel = hl.is_indel(mt.alleles[0], mt.alleles[1]))
    
    return(mt)

def process_consequence_cohort(csq_columns, vcf_metrics_uri):
    ht = hl.import_table(vcf_metrics_uri)

    ht = ht.annotate(locus=hl.parse_variant(ht.ID, reference_genome='GRCh38').locus,
                            alleles=hl.parse_variant(ht.ID, reference_genome='GRCh38').alleles)

    numeric = ['QUAL', 'LEN', 'AC', 'AF', 'AN', 'BaseQRankSum',
           'DP', 'FS', 'MLEAC', 'MLEAF', 'MQ', 'MQRankSum',
           'QD', 'ReadPosRankSum', 'SOR', 'VQSLOD', 'cohort_AC',
           'AB_sample', 'DP_sample', 'GQ_sample', 'VAF_sample',
           'AB_father', 'DP_father', 'GQ_father', 'VAF_father',
           'AB_mother', 'DP_mother', 'GQ_mother', 'VAF_mother', 'GQ_mean']

    mt = ht.to_matrix_table_row_major(columns=numeric, entry_field_name='metrics', col_field_name='SAMPLE')
    mt = mt.key_rows_by(mt.locus, mt.alleles)

    transcript_consequences = mt.CSQ.replace("\[",'').replace("\]",'').replace("\'",'').replace(' ','').split(',').map(lambda x: x.split('\|'))

    transcript_consequences_strs = transcript_consequences.map(lambda x: hl.if_else(hl.len(x)>1, hl.struct(**
                                                           {col: x[i] if col!='Consequence' else x[i].split('&')  
                                                            for i, col in enumerate(csq_columns)}), 
                                                            hl.struct(**{col: '.' if col!='Consequence' else hl.array(['.'])  
                                                            for i, col in enumerate(csq_columns)})))

    mt=mt.annotate_rows(vep=hl.Struct(transcript_consequences = transcript_consequences_strs))

    mt_filt = filter_vep_to_canonical_transcripts(mt)

    mt_csq = process_consequences(mt_filt)

    mt_csq_annotated = add_vep_annotations(mt_csq)

    mt_no_vep = mt_csq_annotated.drop('vep')

    vcf_metrics = mt_no_vep.make_table().to_pandas()

    vcf_metrics.columns = vcf_metrics.columns.str.split('.').str[0]

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

df = process_consequence_cohort(csq_columns, vcf_metrics_uri)
df['isCoding'] = df.Consequence.isin(coding_variants)
df.to_csv(f"{os.path.basename(vcf_metrics_uri).split('.tsv')[0]}_prioritized_csq.tsv", sep='\t',index=False)