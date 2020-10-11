import sys
import pysam
import truvari

in_vcf = pysam.VariantFile(sys.argv[1])
header = in_vcf.header
out_vcf = pysam.VariantFile("/dev/stdout", 'w', header=header)

for entry in in_vcf:
    entry = truvari.copy_entry(entry, header)
    entry.ref = "ATCGATACT"
    entyr.alts = ["A"]
    out_vcf.write(entry)
