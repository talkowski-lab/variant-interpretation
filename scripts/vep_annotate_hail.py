import hail as hl
import sys

vcf_file = sys.argv[1]
vep_annotated_vcf_name = sys.argv[2]

mt = hl.import_vcf(vcf_file, force_bgz=True, reference_genome='GRCh38')
mt = hl.vep(mt, config='vep_config.json', csq=True)
hl.export_vcf(mt, vep_annotated_vcf_name)