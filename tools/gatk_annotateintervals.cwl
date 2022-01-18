cwlVersion: v1.1
class: CommandLineTool
id: gatk_annotateintervals
doc: >-
  Annotates intervals with GC content, and optionally, mappability and segmental-duplication content.
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'broadinstitute/gatk:4.2.4.1'
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      $( inputs.do_explicit_gc_correction ? 'gatk' : 'echo gatk' )
  - position: 1
    shellQuote: false
    valueFrom: >-
      AnnotateIntervals
  - position: 2
    shellQuote: false
    prefix: "-O"
    valueFrom: >-
      ${var prefix = inputs.output_prefix ? inputs.output_prefix : 'intervals'; return prefix + '.annotated.tsv';}
inputs:
  do_explicit_gc_correction:
    type: 'boolean'
    doc: "Choose whether to execute this app. This argument is GATK CNV Best Practice requirement."
  reference:
    type: 'File'
    doc: "Reference fasta"
    secondaryFiles: ['.fai']
    inputBinding:
      position: 2
      prefix: "-R"
  sequence_dictionary:
    type: 'File?'
    doc: "Use the given sequence dictionary as the master/canonical sequence dictionary. Must be a .dict file."
    inputBinding:
      position: 2
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
      position: 2
  input_interval_list:
    type: 'File?'
    secondaryFiles: [{ pattern: ".tbi", required: false }]
    doc: "One or more genomic intervals over which to operate. Use this input when providing interval list files or other file based inputs."
    inputBinding:
      position: 2
      prefix: "-L"
  interval_padding:
    type: 'int?'
    doc: "Amount of padding (in bp) to add to each interval you are including."
    inputBinding:
      position: 2
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
      position: 2
  input_exclude_interval_list:
    type: 'File?'
    secondaryFiles: [{ pattern: ".tbi", required: false }]
    doc: "One or more genomic intervals to exclude from processing. Use this input when providing interval list files or other file based inputs."
    inputBinding:
      position: 2
      prefix: "-XL"
  interval_exclusion_padding:
    type: 'int?'
    doc: "Amount of padding (in bp) to add to each interval you are excluding."
    inputBinding:
      position: 2
      prefix: "-ixp"
  interval_merging_rule:
    type:
      - 'null'
      - type: enum
        name: interval_merging_rule
        symbols: ["ALL","OVERLAPPING_ONLY"]
    doc: "By default, the program merges abutting intervals (i.e. intervals that are directly side-by-side but do not actually overlap) into a single continuous interval. However you can change this behavior if you want them to be treated as separate intervals instead."
    inputBinding:
      position: 2
      prefix: "-imr"
  mappability_track:
    type: 'File?'
    doc: "Path to Umap single-read mappability track in .bed or .bed.gz format (see https://bismap.hoffmanlab.org/). Overlapping intervals must be merged."
    secondaryFiles: ['.tbi']
    inputBinding:
      position: 2
      prefix: "--mappability-track"
  segmental_duplication_track:
    type: 'File?'
    doc: "Path to segmental-duplication track in .bed or .bed.gz format. Overlapping intervals must be merged."
    secondaryFiles: ['.tbi']
    inputBinding:
      position: 2
      prefix: "--segmental-duplication-track"
  feature_query_lookahead:
    type: 'int?'
    doc: "Number of bases to cache when querying feature tracks"
    inputBinding:
      position: 2
      prefix: "--feature-query-lookahead"
  output_prefix:
    type: 'string?'
    doc: "String to use as the prefix for the outputs."
outputs:
  output: { type: 'File?', outputBinding: { glob: '*.annotated.tsv' } }
