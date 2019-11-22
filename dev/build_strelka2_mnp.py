import sys
import gzip
import pdb

mnps_txt = open(sys.argv[1])
strelka2_fn = sys.argv[2]
strelka2_snps = gzip.open(strelka2_fn)

# store mnps, with chrom\tpos as key, list of mnps as values
mnp_dict = {}
head = next(mnps_txt)
for line in mnps_txt:
    info = line.rstrip("\n").split("\t")
    key = info[0] + "\t" + info[1]
    if key not in mnp_dict:
        mnp_dict[key] = []
    mnp_dict[key].append(info[2] + "\t" + info[3])
mnps_txt.close()

header = []
# store header
for line in strelka2_snps:
    header.append(line.decode())
    if line.decode()[0:6] == "#CHROM":
        break
# on rescan, start at data skipping over header
last_pos = strelka2_snps.tell()
# use chrom\tpos as starting point to eval a possible overlapping mnp
new_mnps = {}

for line in strelka2_snps:
    info = line.decode().rstrip("\n").split("\t")
    dkey = info[0] + "\t" + info[1]
    if dkey in mnp_dict:
        max_len = len(max(mnp_dict[dkey], key=len))
        candidates = {}
        # if possible, build longest supported mnp
        ref = info[3]
        alt = info[4]
        cur = ref + "\t" + alt
        pos = info[1]
        for val in mnp_dict[dkey]:
            candidates[val] = 0
        try:
            f = 0
            c = 0
            while f == 0:
                # keep track of current search position in file, if condition met where the next position isn't continous, or exceeds the longest seen mnp, pick up from there in the loop
                cur_pos = strelka2_snps.tell()
                new_in = next(strelka2_snps)
                new_info = new_in.decode().rstrip("\n").split("\t")
                if new_info[0] == info[0] and int(new_info[1]) == (int(pos) + 1):
                    # pdb.set_trace()
                    ref += new_info[3]
                    alt += new_info[4]
                    cur = ref + "\t" + alt
                    pos = new_info[1]
                    if cur in candidates:
                        candidates[cur] = 1
                        c = 1
                    if len(cur) == max_len:
                        f = 1
                        strelka2_snps.seek(cur_pos)
                else:
                    f = 1
                    strelka2_snps.seek(cur_pos)
            if c == 1:
                new_mnps[dkey] = cur
                sys.stderr.write("Candidate mnp found: " + dkey + "\t" + cur + "\n")
                    
        except Exception as e:
            sys.stderr.write(str(e) + "\n")
            sys.stderr.write("Probably hit EOF\n")
            if c == 1:
                new_mnps[dkey] = cur
strelka2_snps.close()
strelka2_modded = open("strelka2_snp_mnp.vcf", "w")
strelka2_modded.write("".join(header))
strelka2_snps = gzip.open(strelka2_fn)
strelka2_snps.seek(last_pos)

# re-scan strelka2 snps, if position in candidate mnps, use that info, else output snp info
for line in strelka2_snps:
    info = line.decode().rstrip("\n").split("\t")
    dkey = info[0] + "\t" + info[1]
    if dkey in new_mnps:
        (ref, alt) = new_mnps[dkey].split("\t")
        strelka2_modded.write(dkey + "\t" + info[2] + "\t" + new_mnps[dkey] + "\t" + "\t".join(info[5:]) + "\n")
        for i in range(1, len(alt), 1):
            skip = next(strelka2_snps)
    else:
        strelka2_modded.write(line.decode())
strelka2_modded.close()