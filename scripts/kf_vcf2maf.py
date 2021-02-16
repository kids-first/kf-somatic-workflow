""" 
Convert vcf to maf
"""

import argparse
import pysam
import pdb

parser = argparse.ArgumentParser()
parser.add_argument('--input_vcf', help='VCF file to subset')
parser.add_argument('--csq_field', help='Key used to find consenquences, usually CSQ or ANN')
parser.add_argument('--maf_vcf_dict', help='MAF descriptor tsv table with MAF Field, Source, VCF Field, VCF Key, VCF Subkey')
parser.add_argument('--tumor_id', help='Sample ID of tumor/affected sample')
parser.add_argument('--normal_id', help='Sample ID of normal/control sample')
parser.add_argument('--output_basename', help='Output maf file base')
args = parser.parse_args()


def get_vcf_field(dargs):
    try:
        maf_field, field, key, subkey = dargs
        if field in id_to_pysam:
            field = id_to_pysam[field]
        vcf_field = getattr(record, field)
        # Most categorized info will come from PICKed consenquence
        if key == "CSQ":
            for csq in vcf_field[key]:
                annot_info = csq.split("|")
                if annot_info[p_idx] == "1":
                    return annot_info[csq_fields.index(subkey)]

        elif key != "":
            vcf_field = vcf_field[key]
            if subkey != "":
                return str(vcf_field[subkey])
        return str(vcf_field)
    except Exception as e:
        print(e)
        return "NA"
        # hold = 1


def slap_on_meta(dargs):
    return "placeholder"


def calc_field(dargs):
    return "placeholder"


in_vcf = pysam.VariantFile(args.input_vcf, threads=8)
# open table and mark which fields are present in vcf to avoid errors
mv_file = open(args.maf_vcf_dict)
head = next(mv_file)
header = head.rstrip('\n').split('\t')
maf_head = []
maf_action = {}
# parse file and assign an action based on where the data would come from
for line in mv_file:
    (maf_field, src, vcf_field, vcf_key, vcf_subkey) = line.rstrip('\n').split('\t')
    maf_head.append(maf_field)
    maf_action[maf_field] = {}
    if src == "VCF Field":
        maf_action[maf_field] = ['get_vcf_field', maf_field, vcf_field, vcf_key, vcf_subkey]
    elif src == "Metadata":
        maf_action[maf_field] = ['slap_on_meta', maf_field]
    else:
        maf_action[maf_field] = ['calc_field', maf_field]

sample_list = list(in_vcf.header.samples)
out_fn = args.output_basename + ".maf"
out_write = open(out_fn, "w")
out_write.write("\t".join(maf_head) + "\n")
# Convert description ID to pysam object reference name
id_to_pysam = {"INFO": "info", "FORMAT": "formats", "FILTER": "filter"}
csq_fields = in_vcf.header.info['CSQ'].description.replace("Consequence annotations from Ensembl VEP. Format: ", "").split("|")
# get index of PICK field
p_idx = csq_fields.index("PICK")

for record in in_vcf:
    out_str = []
    for mfield in maf_head:
        func, dargs = maf_action[mfield][0], maf_action[mfield][1:]
        out_str.append(locals()[func](dargs))
    out_write.write("\t".join(out_str) + "\n")
out_write.close()