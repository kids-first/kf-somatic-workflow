import pysam
import sys
from math import sqrt
import pdb


def get_counts(region, pos, ref_in, alt_in):
    # require min base qual 30, read is paired, map qual > 0
    #MQ formula is sqrt(sum(x_n^2)/n)
    
    mq_vals = []
    dp = 0
    ref_ct = 0
    alt_ct = 0
    for read in region:
        # pdb.set_trace()
        alt_idx = read.positions.index(pos)
        if read.mapq > min_mq \
        and read.query_alignment_qualities[alt_idx] >= min_bq \
        and read.is_proper_pair and not read.is_secondary:
            dp += 1
            mq_vals.append((read.mapq * read.mapq))
            if read.query_alignment_sequence[alt_idx] == ref_in:
                ref_ct += 1
            elif read.query_alignment_sequence[alt_idx] == alt_in:
                alt_ct += 1
    mq = sqrt(sum(mq_vals)/len(mq_vals))
    pdb.set_trace()
    hold = 1
    return []


def get_read_info(contig, pos, ref, alt):
    tum_region = tum_cram_in.fetch(contig, pos, pos)
    tum_ad, tum_dp, tum_gt, tum_mq = get_counts(tum_region, pos, ref, alt)
    norm_region = norm_cram_in.fetch(contig, pos, pos)
    return []

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
    norm_ad, norm_dp, morm_gt, norm_mq) = get_read_info(record.contig, record.pos, record.ref, record.alts[0])