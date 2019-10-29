cwlVersion: v1.0
class: CommandLineTool
id: python_vardict_interval_split
label: Create intervals for VarDict
doc: "This tool takes in an interval list with the WGS coords split by N regions.
It splits the bed files into bed files with a max total number of base bairs, unless the regions is already larger, it stays in it's own file.
Then within the split lists, intervals are split the specified chunks for easier processing by vardict.  This method prevents FP calls caused by
regions with valid ACGT bases from being split between interval lists.  For example, for hg38 canonical chromosomes, using bp_target=60000000 and
intvl_target_size=20000 will yield about 55 bed files, each with about 60M bp worth of coverage (unless the interval was already larger, it will be in it's own list),
split into 20kb chunks."
requirements:
  - class: DockerRequirement
    dockerPull: 'kfdrc/python:2.7.13'
  - class: InlineJavascriptRequirement
baseCommand: [python, -c]
arguments:
  - position: 0
    valueFrom: >-
        def main():
            import sys
            bp_target = $(inputs.bp_target)
            intvl_target_size = $(inputs.intvl_target_size)
            bed_file = open("$(inputs.wgs_bed_file.path)")

            i=0
            intvl_set = {}
            cur_size = 0
            for cur_intvl in bed_file:
                f = 0
                if i not in intvl_set:
                    intvl_set[i] = []
                data = cur_intvl.rstrip("\n").split("\t")
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
            bed_file.close()

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

        if __name__ == "__main__":
            main()

inputs:
  wgs_bed_file: {type: File, doc: "Should be a bed file of WGS regions with N's removed.  GATK calling regions is a good source."}
  bp_target: {type: ['null', int], doc: "Intended max number of base pairs per file.  Existing intervals large than this will NOT be split into another file.", default: 60000000}
  intvl_target_size: {type: ['null', int], doc: "For each file, split each interval into chuck of this size", default: 20000}

outputs:
  split_intervals_bed:
    type: 'File[]'
    outputBinding:
      glob: '*.bed'
