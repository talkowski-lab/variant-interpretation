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

def filter_mt(mt, filter_csq=True, filter_impact=True):
    '''
    mt: can be trio matrix (tm) or matrix table (mt) but must be transcript-level, not variant-level
    '''
    # filter by Consequence
    if filter_csq:
        exclude_csqs = ['intergenic_variant', 'upstream_gene_variant', 'downstream_gene_variant',
                        'synonymous_variant', 'coding_sequence_variant', 'sequence_variant']
        mt = mt.filter_rows(hl.set(exclude_csqs).intersection(
            hl.set(mt.vep.transcript_consequences.Consequence)).size()!=hl.set(mt.vep.transcript_consequences.Consequence).size())

    # filter only canonical transcript
    mt = mt.filter_rows(mt.vep.transcript_consequences.CANONICAL=='YES')

    # filter by Impact and splice/noncoding consequence
    if filter_impact:
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

gnomad_fields = [x for x in list(mt.vep.transcript_consequences[0]) if 'gnomAD' in x]
mt = mt.annotate_rows(all_csqs=hl.set(hl.flatmap(lambda x: x, mt.vep.transcript_consequences.Consequence)),  
                             gnomad_popmax_af=hl.max([hl.or_missing(hl.array(hl.set(mt.vep.transcript_consequences[gnomad_field]))[0]!='',
                                    hl.float(hl.array(hl.set(mt.vep.transcript_consequences[gnomad_field]))[0])) 
                             for gnomad_field in gnomad_fields]))

# Phasing
tmp_ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
tmp_ped.to_csv(f"{prefix}.ped", sep='\t', index=False)
pedigree = hl.Pedigree.read(f"{prefix}.ped", delimiter='\t')

tm = hl.trio_matrix(mt, pedigree, complete_trios=False)
phased_tm = hl.experimental.phase_trio_matrix_by_transmission(tm, call_field='GT', phased_call_field='PBT_GT')

# Mendel errors
all_errors, per_fam, per_sample, per_variant = hl.mendel_errors(mt['GT'], pedigree)
phased_tm = phased_tm.annotate_rows(mendel_code=all_errors.key_by('locus','alleles')[phased_tm.row_key].mendel_code)

# Output 1: grab ClinVar only
clinvar_tm = phased_tm.filter_rows((phased_tm.info.CLNSIG[0].matches('Pathogenic') | phased_tm.info.CLNSIG[0].matches('pathogenic')))
clinvar_tm = clinvar_tm.filter_entries((clinvar_tm.proband_entry.GT.is_non_ref()) | 
                                   (clinvar_tm.mother_entry.GT.is_non_ref()) |
                                   (clinvar_tm.father_entry.GT.is_non_ref()))
clinvar_tm = clinvar_tm.annotate_rows(variant_type='ClinVar_P/LP')
clinvar_tm = clinvar_tm.explode_rows(clinvar_tm.vep.transcript_consequences)
clinvar_tm = filter_mt(clinvar_tm, filter_csq=False, filter_impact=False)
clinvar_tm = get_transmission(clinvar_tm)

# filter out ClinVar benign
mt = mt.filter_rows((hl.is_missing(mt.info.CLNSIG)) |
    ~(mt.info.CLNSIG[0].matches('Benign') | mt.info.CLNSIG[0].matches('benign')))

# filter PASS
mt = mt.filter_rows(mt.filters.size()==0)

# filter out variants containing only these consequences
exclude_csqs = ['intergenic_variant', 'upstream_gene_variant', 'downstream_gene_variant',
                'synonymous_variant', 'coding_sequence_variant', 'sequence_variant']

mt = mt.filter_rows(hl.set(exclude_csqs).intersection(mt.all_csqs).size()!=mt.all_csqs.size())

# filter by AC and gnomAD AF
mt = mt.filter_rows(mt.info.cohort_AC<=ac_threshold)
mt = mt.filter_rows((mt.gnomad_popmax_af<=gnomad_af_threshold) | (hl.is_missing(mt.gnomad_popmax_af)))

# export intermediate VCF
hl.export_vcf(mt, prefix+'_clinical.vcf.bgz', metadata=header)

# export ClinVar TSV
clinvar_tm.entries().flatten().export(prefix+'_clinvar_variants.tsv.gz', delimiter='\t')
