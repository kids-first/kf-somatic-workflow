import pysam
import sys
from math import sqrt
import pdb


def process_snp(snp_read, ct_list, alt_idx, rlen, alen, ref_in, alt_in):
    if snp_read.query_alignment_sequence[alt_idx] == ref_in:
        ct_list[0] += 1
    elif snp_read.query_alignment_sequence[alt_idx] == alt_in:
        ct_list[1] += 1
    return ct_list

def process_mnp(mnp_read, ct_list, alt_idx, rlen, alen, ref_in, alt_in):
    start = alt_idx
    end = alt_idx + rlen
    if mnp_read.query_alignment_sequence[start:end] == ref_in:
        ct_list[0] += 1
    elif mnp_read.query_alignment_sequence[start:end] == alt_in:
        ct_list[1] += 1
    return ct_list

def process_indel(indel_read_read, ct_list, alt_idx, rlen, alen, ref_in, alt_in):
    pdb.set_trace()
    hold=1
    return ct_list
def get_counts(region, vcf_start, vcf_stop, ref_in, alt_in):
    # require read is paired, map qual > 0
    #MQ formula is sqrt(sum(x_n^2)/n)
    
    mq_vals = []
    dp = 0
    ref_alt_ct = [0, 0]
    alt_ct = 0
    call_type = "snp"
    alen = len(alt_in)
    rlen = len(ref_in)
    if alen > 1 or rlen > 1:
        if alen == rlen:
            call_type = "mnp"
    else:
        call_type = "indel"

    for read in region:
        # pysam ref positions are 0-based, need to sutract 1 from vcf start to make compatible
        try:
            alt_idx = read.get_reference_positions().index(vcf_start)
        except ValueError:
            sys.stderr.write('Position not covered by read ' + read.query_name + ' at vcf position ' + read.reference_name + ' ' + str(vcf_start+1) + '\n')
        # get lengths of ref and alts to ensure while seq match
        if read.mapq > min_mq \
        and read.is_proper_pair and not read.is_secondary:
            dp += 1
            ref_alt_ct = globals()["process_" + call_type](read, ref_alt_ct, alt_idx, rlen, alen, ref_in, alt_in)
            mq_vals.append((read.mapq * read.mapq))
    mq = sqrt(sum(mq_vals)/len(mq_vals))
    pdb.set_trace()
    hold = 1
    return ref_alt_ct, dp, mq


def get_read_info(vcf_record):
    contig, start, stop, ref, alt = vcf_record.contig, vcf_record.start, vcf_record.stop, vcf_record.ref, vcf_record.alts[0]
    tum_region = tum_cram_in.fetch(contig, start, stop)
    tum_ad, tum_dp, tum_mq = get_counts(tum_region, start, stop, ref, alt)
    norm_region = norm_cram_in.fetch(contig, start, stop)
    return []


if(len(sys.argv) == 1):
    sys.stderr.write("Needs {fasta reference} {tumor cram} {normal cram} {vcf}\n")
    exit(1)
ref_fasta_fn = sys.argv[1]
tum_cram_fn = sys.argv[2]
norm_cram_fn = sys.argv[3]
vcf_fn = sys.argv[4]

min_bq = 30
min_mq = 0

tum_cram_in = pysam.AlignmentFile(tum_cram_fn, "rc", threads=8, reference_filename=ref_fasta_fn)
norm_cram_in = pysam.AlignmentFile(tum_cram_fn, "rc", threads=8, reference_filename=ref_fasta_fn)
vcf_in = pysam.VariantFile(vcf_fn, threads=8)

for record in vcf_in:
    (ad_i, af_i, dp_i,
    tum_ad, tum_dp, tum_gt, tum_mq,
    norm_ad, norm_dp, morm_gt, norm_mq) = get_read_info(record)