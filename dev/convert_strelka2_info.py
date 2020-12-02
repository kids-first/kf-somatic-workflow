import pysam
import sys
from pysam import VariantFile
import pdb


def _tumor_normal_genotypes(record):
    """Retrieve standard 0/0, 0/1, 1/1 style genotypes from INFO field.
    Normal -- NT field (ref, het, hom, conflict)
    Tumor -- SGT field
      - for SNPs specified as GG->TT for the normal and tumor diploid alleles. These
        can also represent more complex alleles in which case we set at heterozygotes
        pending longer term inclusion of genotypes in Strelka2 directly
        (https://github.com/Illumina/strelka/issues/16)
      - For indels, uses the ref, het, hom convention
    ref: The REF allele from a VCF line
    alt: A list of potentially multiple ALT alleles (rec.ALT.split(";"))
    info: The VCF INFO field
    fname, coords: not currently used, for debugging purposes
    """
    known_names = set(["het", "hom", "ref", "conflict"])
    def name_to_gt(val):
        if val.lower() == "het":
            return "0/1"
        elif val.lower() == "hom":
            return "1/1"
        elif val.lower() in set(["ref", "conflict"]):
            return "0/0"
        else:
            # Non-standard representations, het is our best imperfect representation
            print(record.pos, record.ref, record.alts, val)
            sys.stderr.write("WARN: Non standard representation found, set as 0/1\t" + "\t".join([record.contig, str(record.pos), record.info['NT'], record.info['SGT']]) + "\n")
            return "0/1"
    def alleles_to_gt(val):
        gt_indices = {gt.upper(): i for i, gt in enumerate([record.ref] + list(record.alts))}
        tumor_gts = [gt_indices[x.upper()] for x in val if x in gt_indices]
        if tumor_gts and val not in known_names:
            if max(tumor_gts) == 0:
                tumor_gt = "0/0"
            elif 0 in tumor_gts:
                tumor_gt = "0/%s" % min([x for x in tumor_gts if x > 0])
            else:
                tumor_gt = "%s/%s" % (min(tumor_gts), max(tumor_gts))
        else:
            tumor_gt = name_to_gt(val)
        return tumor_gt
    nt_val = record.info['NT']
    normal_gt = name_to_gt(nt_val)
    sgt_val = record.info['SGT'].split("->")[-1]
    if not sgt_val:
        tumor_gt = "0/0"
    else:
        sgt_val = sgt_val[0].split("->")[-1]
        tumor_gt = alleles_to_gt(sgt_val)
    return tumor_gt, normal_gt
if len(sys.argv) == 1:
    sys.stderr.write('Need input vcf\n')
    exit()
strelka2_in = VariantFile(sys.argv[1])
# create new header adding: ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
strelka2_in.header.add_meta('FORMAT', items=[('ID','GT'), ('Number',1), ('Type','String'), ('Description', 'Genotype')])
updated_vcf = VariantFile('-', 'w', header=strelka2_in.header)
norm_idx = 0
tum_idx = 1
samp_list = list(strelka2_in.header.samples)

for rec in strelka2_in.fetch():
    tumor_gt, normal_gt = _tumor_normal_genotypes(rec)
    pdb.set_trace()
    rec.samples[samp_list[norm_idx]]['GT'] = normal_gt
    rec.samples[samp_list[tum_idx]]['GT'] = tumor_gt
    updated_vcf.write(rec)
updated_vcf.close()


