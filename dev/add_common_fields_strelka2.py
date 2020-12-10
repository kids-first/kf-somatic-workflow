import pysam
import sys
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
        tumor_gt = alleles_to_gt(sgt_val)
    # convert to tuple for use in pysam
    tumor_gt = tuple(map(int, tumor_gt.split("/")))
    normal_gt = tuple(map(int,normal_gt.split ("/")))
    return tumor_gt, normal_gt

def _calc_AD(record):
    """From strelka2 docs:
    SNPs:
    refCounts = Value of FORMAT column $REF + “U” (e.g. if REF="A" then use the value in FORMAT/AU)
    altCounts = Value of FORMAT column $ALT + “U” (e.g. if ALT="T" then use the value in FORMAT/TU)
    tier1RefCounts = First comma-delimited value from $refCounts
    tier1AltCounts = First comma-delimited value from $altCounts
    Indels:
    tier1RefCounts = First comma-delimited value from FORMAT/TAR
    tier1AltCounts = First comma-delimited value from FORMAT/TIR
    """
    tum_alt_ct = 0
    tum_ref_ct = 0
    norm_alt_ct = 0
    norm_ref_ct = 0
    # only indels should have TIR field
    if 'TIR' in record.samples[samp_list[tum_idx]]:
        tum_ref_ct = record.samples[samp_list[tum_idx]]['TAR'][0]
        tum_alt_ct = record.samples[samp_list[tum_idx]]['TIR'][0]
        norm_ref_ct = record.samples[samp_list[norm_idx]]['TAR'][0]
        norm_alt_ct = record.samples[samp_list[norm_idx]]['TIR'][0]

    else:
        snp_ref_key = record.ref + "U"
        snp_alt_key = record.alts[0] + "U"
        tum_ref_ct = record.samples[samp_list[tum_idx]][snp_ref_key][0]
        tum_alt_ct = record.samples[samp_list[tum_idx]][snp_alt_key][0]
        norm_ref_ct = record.samples[samp_list[norm_idx]][snp_ref_key][0]
        norm_alt_ct = record.samples[samp_list[norm_idx]][snp_alt_key][0]
    tum_ad = (tum_ref_ct, tum_alt_ct)
    norm_ad = (norm_ref_ct, norm_alt_ct)
    # return as tuple for pysam compatibility
    return tum_ad, norm_ad



def _calc_AF_AD(record):
    """From strelka2 docs:
    SNPs:
    Somatic allele freqeuncy is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)
    Indels:
    Somatic allele freqeuncy is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)
    """
    tum_af = 0.0
    norm_af = 0.0
    tum_ad, norm_ad = _calc_AD(record)
    try:
        tum_af = tum_ad[1]/sum(list(tum_ad))
    except ZeroDivisionError:
        sys.stderr.write(str(e) + "0 div error while calculating af for tumor at " + record.contig + " " 
        + str(record.pos) + ", setting to 0.0\n")
    except Exception as e:
        sys.stderr.write(str(e) + "\nError while calculating af for tumor at " + record.contig + " " 
        + str(record.pos) + "\n")
    try:
        norm_af = norm_ad[1]/sum(list(norm_ad))
    except ZeroDivisionError:
        sys.stderr.write(str(e) + "0 div error while calculating af for normal at " + record.contig + " " 
        + str(record.pos) + ", setting to 0.0\n")
    except Exception as e:
        sys.stderr.write(str(e) + "\nError while calculating af for normal at " + record.contig + " " 
        + str(record.pos) + "\n")

    return tum_af, norm_af, tum_ad, norm_ad


if len(sys.argv) == 1:
    sys.stderr.write('Need input vcf and output file prefix\n')
    exit()
strelka2_in = pysam.VariantFile(sys.argv[1])
out_fn = sys.argv[2] + ".vcf.gz"
# create new header adding: ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
# and: ##FORMAT=<ID=AF,Number=.,Type=Float,Description="Allele frequency, as recommended by strelka2 docs: <ALT>U/<REF>U+<ALT>U (somatic snps), 'TIR/TIR+TAR (somatic indels)"
# and: ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed. Added in post for compatibility"
strelka2_in.header.add_meta('FORMAT', items=[('ID','GT'), ('Number',1), ('Type','String'), ('Description', 'Genotype converted for cross-compatibility using bcbio method')])
strelka2_in.header.add_meta('FORMAT', items=[('ID','AF'), ('Number', '.'), ('Type', 'Float'), ('Description', 'Allele frequency, as recommended by strelka2 docs: <ALT>U/<REF>U+<ALT>U (somatic snps), TIR/TIR+TAR (somatic indels)')])
strelka2_in.header.add_meta('FORMAT', items=[('ID','AD'), ('Number', 'R'), ('Type', 'Integer'), ('Description', 'Allelic depths for the ref and alt alleles in the order listed. Added in post for compatibility')])
updated_vcf = pysam.VariantFile(out_fn, 'w', header=strelka2_in.header, threads=8)
norm_idx = 0
tum_idx = 1
samp_list = list(strelka2_in.header.samples)

for rec in strelka2_in.fetch():
    # get new GT, AF, AD values
    tumor_gt, normal_gt = _tumor_normal_genotypes(rec)
    tum_af_val, norm_af_val ,tum_ad_val, norm_ad_val = _calc_AF_AD(rec)
    # Set values in record and print
    rec.samples[samp_list[norm_idx]]['GT'] = normal_gt
    rec.samples[samp_list[tum_idx]]['GT'] = tumor_gt
    rec.samples[samp_list[norm_idx]]['AF'] = norm_af_val
    rec.samples[samp_list[tum_idx]]['AF'] = tum_af_val
    rec.samples[samp_list[norm_idx]]['AD'] = norm_ad_val
    rec.samples[samp_list[tum_idx]]['AD'] = tum_ad_val

    updated_vcf.write(rec)
updated_vcf.close()
pysam.tabix_index(out_fn, preset="vcf")
