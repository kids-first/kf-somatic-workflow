import pysam
import sys
from collections import Counter
import pdb


def get_consensus_gt(sname):
    # get max number of times a GT was seen
    max_val = max(gt_cur[sname].values())
    # see if max was ever seen more than once - meaning a consensus was not reached
    count = Counter(gt_cur[sname].values())

    # get GT values if consensus reached
    gt = ""
    status = "CONFLICT"
    if count[max_val] == 1:
        if len(gt_cur[sname].values()) > 1:
            status = "WARN"
        else:
            status = "UNANIMOUS"
        keys = list(gt_cur[sname].keys())
        vals = list(gt_cur[sname].values())
        gt = keys[vals.index(max_val)]
    return gt, status


if len(sys.argv) == 1:
    sys.stderr.write("Need consensus vcf and vcf csv\n")
    exit(1)

consensus_vcf = pysam.VariantFile(sys.argv[1])

vcf_list = []

for vcf in sys.argv[2].split(','):
    vcf_list.append(pysam.VariantFile(vcf))

print("chrom\tstart\tstop\tref\talt\tsample1 GT\tstatus1\tsample2 GT status2")
for record in consensus_vcf.fetch():
    gt_cur = {}
    s1 = record.samples.keys()[0]
    s2 = record.samples.keys()[1]
    gt_cur[s1] = {}
    gt_cur[s2] = {}
    rec_str = "\t".join([record.contig, str(record.start), str(record.stop), record.ref, record.alts[0]])
    for i in range(len(vcf_list)):
        check = vcf_list[i].fetch(record.contig, record.start, record.stop)
        if check:
            # multi-allelics will pop up displaying a GT adjusted for this. For example, 0/0/1 is likely two hits at the position
            # first is 0/0, second is 0/1. Will track position and adjust what GT to pluck
            j = 0
            for result in check:
                res_str = "\t".join([result.contig, str(result.start), str(result.stop), result.ref, result.alts[0]])
                if rec_str == res_str:
                    gt1 = str(result.samples[s1]['GT'])
                    gt2 = str(result.samples[s2]['GT'])
                    # GT > 2 is a former multi-allelic, shift to capure desired genotype
                    if len(result.samples[s2]['GT']) > 2:
                        gt2 = '(' + str(result.samples[s2]['GT'][(0+j)]) + ', ' + str(result.samples[s2]['GT'][(1+j)]) + ')'

                    if gt1 not in gt_cur[s1]:
                        gt_cur[s1][gt1] = 0
                    if gt2 not in gt_cur[s2]:
                        gt_cur[s2][gt2] = 0
                    gt_cur[s1][gt1] += 1
                    gt_cur[s2][gt2] += 1
                    break
                j += 1
    gt1_val, status1 = get_consensus_gt(s1)
    gt2_val, status2 = get_consensus_gt(s2)
    if status2 == "CONFLICT":
        pdb.set_trace()
        hold = 1
    print(rec_str + "\t" + "\t".join([gt1_val, status1, gt2_val, status2]))
