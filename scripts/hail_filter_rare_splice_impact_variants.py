from pyspark.sql import SparkSession
import hail as hl
import numpy as np
import sys
import ast
import os

vcf_file = sys.argv[1]
output_filename = sys.argv[2]
cores = sys.argv[3]  # string
mem = int(np.floor(float(sys.argv[4])))
ac_threshold = int(sys.argv[5])
gnomad_af_threshold = float(sys.argv[6])

hl.init()

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

# filter out variants containing only these consequences
exclude_csqs = ['intergenic_variant', 'upstream_gene_variant', 'downstream_gene_variant',
                'synonymous_variant', 'coding_sequence_variant', 'sequence_variant']

mt = mt.annotate_rows(all_csqs=hl.set(hl.flatmap(lambda x: x, mt.vep.transcript_consequences.Consequence)),
                             gnomad_af=hl.or_missing(hl.array(hl.set(mt.vep.transcript_consequences.gnomADg_AF))[0]!='', 
                                                  hl.float(hl.array(hl.set(mt.vep.transcript_consequences.gnomADg_AF))[0])))
mt = mt.filter_rows(hl.set(exclude_csqs).intersection(mt.all_csqs).size()!=mt.all_csqs.size())

# filter by AC and gnomAD AF
mt = mt.filter_rows(mt.info.cohort_AC<=ac_threshold)
mt = mt.filter_rows((mt.gnomad_af<=gnomad_af_threshold) | (hl.is_missing(mt.gnomad_af)))

# filter splice variants and MODERATE/HIGH impact variants
keep_vars = ['splice_donor_5th_base_variant', 'splice_region_variant', 'splice_donor_region_variant',  # splice variants
              'non_coding_transcript_exon_variant']

mt = mt.filter_rows(hl.any(lambda csq: hl.array(keep_vars).contains(csq), mt.all_csqs) |
                  (hl.any(lambda impact: hl.array(['HIGH','MODERATE']).contains(impact), 
                          mt.vep.transcript_consequences.IMPACT)))

# filter PASS
mt = mt.filter_rows(mt.filters.size()==0)
hl.export_vcf(mt, output_filename, metadata=header, tabix=True)