import pysam
import sys
from collections import Counter
import pdb


def pysam_gt_to_str(gt_obj, p_flag):
    delim = "/"
    if p_flag:
        delim = "|"
    gt_str = delim.join(map(str, gt_obj))
    return gt_str


def get_consensus_gt(sname):
    # get max number of times a GT was seen
    max_val = max(gt_cur[sname].values())
    # see if max was ever seen more than once - meaning a consensus was not reached
    count = Counter(gt_cur[sname].values())

    # get GT values if consensus reached
    gt = ""
    status = "CONTESTED"
    if count[max_val] == 1:
        if len(gt_cur[sname].values()) > 1:
            status = "MAJORITY"
        else:
            status = "UNANIMOUS"
        keys = list(gt_cur[sname].keys())
        vals = list(gt_cur[sname].values())
        gt = keys[vals.index(max_val)]
    return gt, status


if len(sys.argv) == 1:
    sys.stderr.write("Need consensus vcf, vcf csv, normal ID, tumor ID\n")
    exit(1)
consensus_vcf = pysam.VariantFile(sys.argv[1])
vcf_list = []

for vcf in sys.argv[2].split(','):
    vcf_list.append(pysam.VariantFile(vcf, threads=2))
s1 = sys.argv[3]
s2 = sys.argv[4]
print("chrom\tstart\tstop\tref\talt\tcallers\tnorm GT\tnorm status\ttum original genotypes\ttum GT\ttum status\ttum original genotypes")
m = 1000
rec_ct = 1
for record in consensus_vcf.fetch():
    if rec_ct % m == 0:
        sys.stderr.write("Processed " + str(rec_ct) + " records\n")
        sys.stderr.flush()
    gt_cur = {}
    gt_cur[s1] = {}
    gt_cur[s2] = {}
    rec_str = "\t".join([record.contig, str(record.start), str(record.stop), record.ref, record.alts[0]])
    debug_list = []
    tum_gt_list = []
    norm_gt_list = []
    callers = ",".join(record.info['CALLERS'])
    for i in range(len(vcf_list)):
        check = vcf_list[i].fetch(record.contig, record.start, record.stop)
        if check:
            # multi-allelics will pop up displaying a GT adjusted for this. For example, 0/0/1 is likely two hits at the position
            # first is 0/0, second is 0/1. Will track position and adjust what GT to pluck
            j = 0
            for result in check:
                res_str = "\t".join([result.contig, str(result.start), str(result.stop), result.ref, result.alts[0]])
                if rec_str == res_str:
                    phased = result.samples[s2].phased
                    gt1 = result.samples[s1]['GT']
                    gt2 = result.samples[s2]['GT']
                    gt1_str = pysam_gt_to_str(gt1, phased)
                    norm_gt_list.append(gt1_str)
                    gt2_str = pysam_gt_to_str(gt2, phased)
                    tum_gt_list.append(gt2_str)

                    # GT > 2 is a former multi-allelic, shift to capure desired genotype; likely mutect2
                    if len(gt2) > 2:
                        gt2 = (result.samples[s2]['GT'][(0+j)], result.samples[s2]['GT'][(1+j)])
                        gt2_str = pysam_gt_to_str(gt2, phased)
                    # 0,1 and 1,0 are treated the same ; likely vardict, mutect2 on phased calls
                    if gt2 == (1, 0): # and not phased:
                        gt2_str = '0/1'
                    # For standardizing, will consider 0/1 == 0|1; likely mutect2
                    if phased:
                        gt1_str = gt1_str.replace("|", "/")
                        gt2_str = gt2_str.replace("|", "/")
                    if gt1_str not in gt_cur[s1]:
                        gt_cur[s1][gt1_str] = 0
                    if gt2_str not in gt_cur[s2]:
                        gt_cur[s2][gt2_str] = 0
                    gt_cur[s1][gt1_str] += 1
                    gt_cur[s2][gt2_str] += 1

                    debug_list.append(result)
                    break
                j += 1
    gt1_val, status1 = get_consensus_gt(s1)
    gt2_val, status2 = get_consensus_gt(s2)
    if status2 == "CONTESTED":
        gt2_val = "0/1"
        # sys.stderr.write("DEBUG MODE: OUTPUTTING VARIANTS\n")
        # sys.stderr.write(str(record))
        # for debug_rec in debug_list:
        #     sys.stderr.write(str(debug_rec))

        # pdb.set_trace()
        # hold = 1
    print(rec_str + "\t" + callers + "\t"
    + "\t".join([gt1_val, status1, ",".join(norm_gt_list), gt2_val, status2, ",".join(tum_gt_list)]))
    rec_ct += 1
