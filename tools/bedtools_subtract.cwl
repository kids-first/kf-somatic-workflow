cwlVersion: v1.2
class: CommandLineTool
id: bedtools_subtract
doc: |
  bedtools subtract bedtools subtract searches for features in B that overlap A
  by at least the number of base pairs given by the -f option. If an overlapping
  feature is found in B, the overlapping portion is removed from A and the
  remaining portion of A is reported. If a feature in B overlaps all of a feature
  in A, the A feature will not be reported. If a feature in B does not overlap a
  feature in A by at least the -f amount, the A feature will be reported in its
  entirety.
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.cpu)
  - class: DockerRequirement
    dockerPull: 'staphb/bedtools:2.31.0'
baseCommand: [bedtools, subtract]
stdout: $(inputs.output_filename) 
inputs:
  a_bed: {type: 'File', inputBinding: { prefix: '-a', position: 8 }, doc: "Base bed from which the b_bed will subtract regions." }
  b_bed: {type: 'File', inputBinding: { prefix: '-b', position: 8 }, doc: "BED to subtract from the a_bed." }
  output_filename: { type: 'string?', default: "subtracted.bed", doc: "Name for the output file." }
  remove_any_overlap: { type: 'boolean?', inputBinding: { prefix: '-A', position: 2 }, doc: "Remove entire feature if any overlap.  That is, by default, only subtract the portion of A that overlaps B. Here, if any overlap is found (or -f amount), the entire feature is removed." }
  remove_sum_overlap: { type: 'boolean?', inputBinding: { position: 2, prefix: "-N"}, doc: "Same as -A except when used with -f, the amount is the sum of all features (not any single feature)." }
  strandedness_match: { type: 'boolean?', inputBinding: { position: 2, prefix: "-s"}, doc: "Require same strandedness.  That is, only report hits in B that overlap A on the _same_ strand.  - By default, overlaps are reported without respect to strand." }
  strandedness_mismatch: { type: 'boolean?', inputBinding: { position: 2, prefix: "-S"}, doc: "Require different strandedness.  That is, only report hits in B that overlap A on the _opposite_ strand.  - By default, overlaps are reported without respect to strand." }
  fraction_overlap_a: { type: 'float?', inputBinding: { position: 2, prefix: "-f"}, doc: "Minimum overlap required as a fraction of A.  - Default is 1E-9 (i.e., 1bp).  - FLOAT (e.g. 0.50)" }
  fraction_overlap_b: { type: 'float?', inputBinding: { position: 2, prefix: "-F"}, doc: "Minimum overlap required as a fraction of B.  - Default is 1E-9 (i.e., 1bp).  - FLOAT (e.g. 0.50)" }
  fraction_overlap_a_is_b: { type: 'boolean?', inputBinding: { position: 2, prefix: "-r"}, doc: "Require that the fraction overlap be reciprocal for A AND B.  - In other words, if -f is 0.90 and -r is used, this requires that B overlap 90% of A and A _also_ overlaps 90% of B." }
  either_a_or_b_f: { type: 'boolean?', inputBinding: { position: 2, prefix: "-e"}, doc: "Require that the minimum fraction be satisfied for A OR B.  - In other words, if -e is used with -f 0.90 and -F 0.10 this requires that either 90% of A is covered OR 10% of  B is covered.  Without -e, both fractions would have to be satisfied." }
outputs:
  subtracted_bed:
    type: stdout 
