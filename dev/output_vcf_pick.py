""" 
Quick script to output vcf with PICK only hits
"""

import argparse
import pysam

parser = argparse.ArgumentParser()
parser.add_argument('--input_vcf', help='VCF file to subset')
parser.add_argument('--csq_field', help='Key used to find consenquences, usually CSQ or ANN')
args = parser.parse_args()

in_vcf = pysam.VariantFile(args.input_vcf, threads=8)
csq_fields = in_vcf.header.info['CSQ'].description.replace("Consequence annotations from Ensembl VEP. Format: ", "").split("|")

# get index of PICK field
p_idx = csq_fields.index("PICK")

print("chrom\tpos\tref\talt\tpicked_ct\tpicked_csq_csv")

for record in in_vcf.fetch():
    picked_csqs = []
    for csq in record.info['CSQ']:
        info = csq.split("|")
        if info[p_idx] == '1':
            picked_csqs.append(csq)
    print("\t".join([record.contig, str(record.pos), record.ref, record.alts[0], str(len(picked_csqs)), ",".join(picked_csqs)]))
