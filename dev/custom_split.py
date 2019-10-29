#!/usr/bin/env python3

import sys
# in actual tool, these are inputs
# max number bp per interval file
bp_target = 60000000
# in each interval file, chop up into chunks of this size in bp
intvl_target_size = 20000
bed_file = open(sys.argv[1])

# simple key to hold sets of intervals for output in each file
i = 0
intvl_set = {}
cur_size = 0
for cur_intvl in bed_file:
    f = 0
    if i not in intvl_set:
        intvl_set[i] = []
    data = cur_intvl.rstrip('\n').split('\t')
    (chrom, start, end) = (data[0], data[1], data[2])
    intvl_size = int(end) - int(start)
    # if the interval size is already bigger than the max target, just make it it's own interva list file
    if intvl_size >= bp_target:
        if len(intvl_set[i]) != 0:
            i += 1
            intvl_set[i] = []
            f = 1
    # similar, if adding the interval being processed increases the interval list beyond bp_target, end that list and put the current one in a new list
    elif cur_size + intvl_size > bp_target:
        if len(intvl_set[i]) != 0:
            i += 1
            intvl_set[i] = []
            cur_size = intvl_size
    # if addin g to th ecurrent list is still under target, add to current list
    else:
        cur_size += intvl_size
    intvl_set[i].append([chrom, start, end])
    if f == 1:
        i += 1
        cur_size = 0

for set_i, invtl_list in sorted(intvl_set.items()):
    set_size = 0
    out = open("set_" + str(set_i) + ".bed", "w")
    # for each interval list set in the dict, split intervals in intvl_target_size pieces
    for intervals in invtl_list:
        (chrom, start, end) = (intervals[0], intervals[1], intervals[2])
        intvl_size = int(end) - int(start)
        set_size += intvl_size
        for j in range(int(start), int(end), intvl_target_size):
            new_end = j + intvl_target_size
            if new_end > int(end):
                new_end = end
            out.write(chrom + "\t" + str(j) + "\t" + str(new_end) + "\n")
    # for informational purpose, output total number of bp covered in each list
    sys.stderr.write("Set " + str(set_i) + " size:\t" + str(set_size) + "\n")
    out.close()