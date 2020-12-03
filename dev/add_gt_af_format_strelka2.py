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
        tumor_gt = alleles_to_gt(sgt_val)
    # convert to tuple for use in pysam
    # pdb.set_trace()
    tumor_gt = tuple(map(int, tumor_gt.split("/")))
    normal_gt = tuple(map(int,normal_gt.split ("/")))
    return tumor_gt, normal_gt


def _calc_AF(record):
    """From bcio:
    Strelka2 doesn't report exact AF for a variant, however it can be calculated as alt_counts/dp from existing fields:
    somatic
      snps:    GT:DP:FDP:SDP:SUBDP:AU:CU:GU:TU                 dp=DP                {ALT}U[0] = alt_counts(tier1,tier2)
      indels:  GT:DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50  dp=DP                TIR = alt_counts(tier1,tier2)
    """
    # only indels should have TIR field
    af = 0.0
    if 'TIR' in record.samples[samp_list[tum_idx]]:
        try:
            af = sum(list(record.samples[samp_list[tum_idx]]['TIR']))/record.samples[samp_list[tum_idx]]['DP']
        except Exception as e:
            sys.stderr.write(str(e) + "\nError while calculating af for indel at " + record.contig + " " 
            + str(record.pos))
        return af
    else:
        snp_key = record.alts[0] + "U"
        try:
            af = sum(list(record.samples[samp_list[tum_idx]][snp_key]))/record.samples[samp_list[tum_idx]]['DP']
        except Exception as e:
            sys.stderr.write(str(e) + "\nError while calculating af for snp at " + record.contig + " " 
            + str(record.pos))
        return af




if len(sys.argv) == 1:
    sys.stderr.write('Need input vcf and output file prefix\n')
    exit()
strelka2_in = pysam.VariantFile(sys.argv[1])
out_fn = sys.argv[2] + ".vcf.gz"
# create new header adding: ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
# and: ##FORMAT=<ID=AF,Number=.,Type=Float,Description="Allele frequency, as calculated in bcbio: <ALT>U/DP (somatic snps), 'TIR (somatic indels)"
strelka2_in.header.add_meta('FORMAT', items=[('ID','GT'), ('Number',1), ('Type','String'), ('Description', 'Genotype converted for cross-compatibility using bcbio method')])
strelka2_in.header.add_meta('FORMAT', items=[('ID','AF'), ('Number', '.'), ('Type', 'Float'), ('Description', 'Allele frequency, as calculated in bcbio: <ALT>U/DP (somatic snps), TIR (somatic indels)')])
updated_vcf = pysam.VariantFile(out_fn, 'w', header=strelka2_in.header, threads=8)
norm_idx = 0
tum_idx = 1
samp_list = list(strelka2_in.header.samples)

for rec in strelka2_in.fetch():
    tumor_gt, normal_gt = _tumor_normal_genotypes(rec)
    # some commented values as I decide whether it makes sense to include AF for normal sample
    # tum_af_val, norm_af_val = _calc_AF(rec)
    tum_af_val = _calc_AF(rec)
    rec.samples[samp_list[norm_idx]]['GT'] = normal_gt
    rec.samples[samp_list[tum_idx]]['GT'] = tumor_gt
    # rec.samples[samp_list[norm_idx]]['AF'] = norm_af_val
    rec.samples[samp_list[tum_idx]]['AF'] = tum_af_val
    updated_vcf.write(rec)
updated_vcf.close()
pysam.tabix_index(out_fn, preset="vcf")


