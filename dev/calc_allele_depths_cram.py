import pysam
import sys
from math import sqrt
import pdb


def process_snp(snp_read, ct_list, alt_idx, rlen, alen, ref_in, alt_in):
    if alt_idx is None:
        return ct_list
    if snp_read.query_alignment_sequence[alt_idx] == ref_in:
        ct_list[0] += 1
    elif snp_read.query_alignment_sequence[alt_idx] == alt_in:
        ct_list[1] += 1
    return ct_list


def process_mnp(mnp_read, ct_list, alt_idx, rlen, alen, ref_in, alt_in):
    # pdb.set_trace()
    if alt_idx is None:
        return ct_list
    start = alt_idx
    end = alt_idx + alen
    if mnp_read.query_alignment_sequence[start:end] == alt_in:
        ct_list[1] += 1
    else:
        end = alt_idx + rlen
        if mnp_read.query_alignment_sequence[start:end] == ref_in:
            ct_list[0] += 1
    return ct_list
def process_ins(ins_read, ct_list, alt_idx, rlen, alen, ref_in, alt_in):
    if alt_idx is None:
        return ct_list
    start = alt_idx
    end = alt_idx + alen

    if ins_cigar in ins_read.cigarstring:
        blocks = ins_read.get_blocks()
        for i in range(0, len(blocks)-1, 1):
            # looking for a break that appears in read, end of one seg will have same position of start on the other
            if blocks[i][1] == record.pos:
                if blocks[(i+1)][0] == blocks[i][1] and ins_read.query_alignment_sequence[start:end] == alt_in:
                    ct_list[1] +=1
                    return ct_list
            break
    else:
        end = alt_idx + rlen
        if ins_read.query_alignment_sequence[start:end] == ref_in:
            ct_list[0] += 1
    return ct_list


def process_del(del_read, ct_list, alt_idx, rlen, alen, ref_in, alt_in):
    if del_cigar in del_read.cigarstring:
        # block array of tuples as start and ends of read coverage, if del, diff will equal del _len
        blocks = del_read.get_blocks()
        for i in range(0, len(blocks)-1, 1):
            if blocks[i][1] == record.pos:
                # if true, stop here, if not, could still cover ref, should still check
                if blocks[(i+1)][0] == del_end:
                    ct_list[1] +=1
                    return ct_list
            break
    elif alt_idx:
        start = alt_idx
        end = alt_idx + rlen
        if del_read.query_alignment_sequence[start:end] == ref_in:
            ct_list[0] += 1
    return ct_list


def get_counts(region, vcf_start, vcf_stop, ref_in, alt_in, call_type, rlen, alen):
    # require read is paired, map qual > 0
    #MQ formula is sqrt(sum(x_n^2)/n)

    mq_vals = []
    mq0 = 0
    dp = 0
    ref_alt_ct = [0, 0]

    for read in region:
        mq_vals.append((read.mapq * read.mapq))
        if read.mapq > min_mq \
        and read.is_proper_pair and not read.is_secondary:
            dp += 1
            try:
                alt_idx = read.get_reference_positions().index(vcf_start)
            except ValueError:
                # sys.stderr.write('Position not covered by read ' + read.query_name + ' at vcf position ' + read.reference_name + ' ' + str(vcf_start+1) + " " + ref_in + alt_in + '\n')
                alt_idx = None
            ref_alt_ct = globals()["process_" + call_type](read, ref_alt_ct, alt_idx, rlen, alen, ref_in, alt_in)

        if read.mapq == 0:
            mq0 += 1
    try:
        mq = sqrt(sum(mq_vals)/len(mq_vals))
    except ZeroDivisionError:
        pdb.set_trace()
        hold = 1
    return ref_alt_ct, dp, mq, mq0


def get_read_info(vcf_record):
    contig, start, stop, ref, alt = vcf_record.contig, vcf_record.start, vcf_record.stop, vcf_record.ref, vcf_record.alts[0]
    call_type = "snp"
    alen = len(alt)
    rlen = len(ref)
    if alen > 1 or rlen > 1:
        call_type = "mnp"
        if rlen > 1 and alen == 1:
            call_type = "del"
            # make global so only created once per variant
            global del_cigar
            del_cigar = str(rlen-1) + "D"
            global del_end
            del_end = record.pos + rlen - 1
        elif alen > 1 and rlen == 1:
            call_type = "ins"
            global ins_cigar
            ins_cigar = str(alen-1) + "I"

    tum_region = tum_cram_in.fetch(contig, start, stop)
    tum_ad, tum_dp, tum_mq, tum_mq0 = get_counts(tum_region, start, stop, ref, alt, call_type, rlen, alen)
    norm_region = norm_cram_in.fetch(contig, start, stop)
    norm_ad, norm_dp, norm_mq, norm_mq0 = get_counts(norm_region, start, stop, ref, alt, call_type, rlen, alen)
    return tum_ad, tum_dp, tum_mq, tum_mq0, norm_ad, norm_dp, norm_mq, norm_mq0, call_type


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
norm_cram_in = pysam.AlignmentFile(norm_cram_fn, "rc", threads=8, reference_filename=ref_fasta_fn)
vcf_in = pysam.VariantFile(vcf_fn, threads=8)

processed = 0
m = 100
print("\t".join(["chrom", "position", "ref", "alt", "call type", "tum ref ct", "tum alt ct", "tum depth", "tum mq", "norm ref ct", "norm alt ct", "norm depth", "norm mq"]))
for record in vcf_in:
    if processed % m == 0:
        sys.stderr.write("Processed " + str(processed) + " records\n")
        sys.stderr.flush()
    (tum_ad, tum_dp, tum_mq, tum_mq0, norm_ad, norm_dp, norm_mq, norm_mq0, call_type) = get_read_info(record)
    print_list = [tum_ad[0], tum_ad[1], tum_dp, tum_mq, tum_mq0, norm_ad[0], norm_ad[1], norm_dp, norm_mq, norm_mq0]
    print_list = list(map(str, print_list))
    print(record.contig + "\t" + str(record.pos) + "\t" + record.ref + "\t" + record.alts[0] + "\t" + call_type + "\t" + "\t".join(print_list))
    processed += 1
vcf_in.close()
tum_cram_in.close()
norm_cram_in.close()