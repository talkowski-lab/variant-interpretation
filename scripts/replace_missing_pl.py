from pysam import VariantFile
import sys

vcf_file = sys.argv[1]
out_vcf = sys.argv[2]

bcf_in = VariantFile(vcf_file)  # auto-detect input format
bcf_out = VariantFile(out_vcf, 'w', header=bcf_in.header)

for rec in bcf_in.fetch():
    new_rec = rec.copy()
    for sample in rec.samples:
        if rec.samples[sample]['PL'][0] is None:
            new_rec.samples[sample]['PL'] = (0,0,0)
    bcf_out.write(new_rec)

bcf_in.close()
bcf_out.close()