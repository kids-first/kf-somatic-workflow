cwlVersion: v1.1
class: CommandLineTool
id: gatk_preprocessintervals
doc: "Prepares bins for coverage collection"
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'broadinstitute/gatk:4.2.4.1'
baseCommand: [gatk, PreprocessIntervals]
arguments:
  - position: 1
    shellQuote: false
    prefix: "-O"
    valueFrom: >-
      ${var pre = inputs.output_prefix ? inputs.output_prefix : 'output'; return pre+'.preprocessed_intervals.interval_list'}
inputs:
  output_prefix:
    type: 'string?'
    doc: "String to use as the prefix for the outputs."
  reference:
    type: 'File'
    doc: "Reference fasta"
    secondaryFiles: ['.fai']
    inputBinding:
      prefix: "-R"
  sequence_dictionary:
    type: 'File?'
    doc: "Use the given sequence dictionary as the master/canonical sequence dictionary. Must be a .dict file."
    inputBinding:
      prefix: "--sequence-dictionary"
  input_intervals:
    type:
      - 'null'
      - type: array
        items: string
        inputBinding:
          prefix: '-L'
    doc: "One or more genomic intervals over which to operate. Use this input when providing the intervals as strings."
    inputBinding:
      position: 0
  input_interval_list:
    type: 'File?'
    secondaryFiles: [{ pattern: ".tbi", required: false }]
    doc: "One or more genomic intervals over which to operate. Use this input when providing interval list files or other file based inputs."
    inputBinding:
      prefix: "-L"
  interval_padding:
    type: 'int?'
    doc: "Amount of padding (in bp) to add to each interval you are including."
    inputBinding:
      prefix: "-ip"
  input_exclude_intervals:
    type:
      - 'null'
      - type: array
        items: string
        inputBinding:
          prefix: '-XL'
    doc: "One or more genomic intervals to exclude from processing. Use this input when providing the intervals as strings."
    inputBinding:
      position: 0
  input_exclude_interval_list:
    type: 'File?'
    secondaryFiles: [{ pattern: ".tbi", required: false }]
    doc: "One or more genomic intervals to exclude from processing. Use this input when providing interval list files or other file based inputs."
    inputBinding:
      prefix: "-XL"
  interval_exclusion_padding:
    type: 'int?'
    doc: "Amount of padding (in bp) to add to each interval you are excluding."
    inputBinding:
      prefix: "-ixp"
  padding:
    type: 'int?'
    doc: "Length (in bp) of the padding regions on each side of the intervals."
    inputBinding:
      prefix: "--padding"
  bin_length:
    type: 'int?'
    inputBinding:
      prefix: "--bin-length"
  interval_merging_rule:
    type:
      - 'null'
      - type: enum
        name: interval_merging_rule
        symbols: ["ALL","OVERLAPPING_ONLY"]
    doc: "By default, the program merges abutting intervals (i.e. intervals that are directly side-by-side but do not actually overlap) into a single continuous interval. However you can change this behavior if you want them to be treated as separate intervals instead."
    inputBinding:
      prefix: "-imr"
  interval_set_rule:
    type:
      - 'null'
      - type: enum
        name: interval_set_rule
        symbols: ["UNION","INTERSECTION"]
    doc: "By default, the program will take the UNION of all intervals specified using -L and/or -XL. However, you can change this setting for -L, for example if you want to take the INTERSECTION of the sets instead. E.g. to perform the analysis only on chromosome 1 exomes, you could specify -L exomes.intervals -L 1 --interval-set-rule INTERSECTION. However, it is not possible to modify the merging approach for intervals passed using -XL (they will always be merged using UNION). Note that if you specify both -L and -XL, the -XL interval set will be subtracted from the -L interval set."
    inputBinding:
      prefix: "-isr"
outputs:
  output: { type: 'File', outputBinding: { glob: '*interval_list' } }
