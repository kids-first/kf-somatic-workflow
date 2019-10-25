#!/usr/bin/env python3

import sys

bp_target = 60000000
intvl_target_size = 20000
bed_file = open(sys.argv[1])

i=0
intvl_set = {}
cur_size = 0
for cur_intvl in bed_file:
    f = 0
    if i not in intvl_set:
        intvl_set[i] = []
    data = cur_intvl.rstrip('\n').split('\t')
    (chrom, start, end) = (data[0], data[1], data[2])
    intvl_size = int(end) - int(start)
    if intvl_size >= bp_target:
        if len(intvl_set[i]) != 0:
            i += 1
            intvl_set[i] = []
            f = 1
    elif cur_size + intvl_size > bp_target:
        if len(intvl_set[i]) != 0:
            i += 1
            intvl_set[i] = []
            cur_size = intvl_size
    else:
        cur_size += intvl_size
    intvl_set[i].append([chrom, start, end])
    if f == 1:
        i += 1
        cur_size = 0

for set_i, invtl_list in sorted(intvl_set.items()):
    set_size = 0
    out = open("set_" + str(set_i) + ".bed", "w")
    for intervals in invtl_list:
        (chrom, start, end) = (intervals[0], intervals[1], intervals[2])
        intvl_size = int(end) - int(start)
        set_size += intvl_size
        for j in range(int(start), int(end), intvl_target_size):
            new_end = j + intvl_target_size
            if new_end > int(end):
                new_end = end
            out.write(chrom + "\t" + str(j) + "\t" + str(new_end) + "\n")
    sys.stderr.write("Set " + str(set_i) + " size:\t" + str(set_size) + "\n")
    out.close()