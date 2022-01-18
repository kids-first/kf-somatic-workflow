cwlVersion: v1.1
class: CommandLineTool
id: gatk_collectreadcounts
doc: "Collects read counts at specified intervals. The count for each interval is calculated by counting the number of read starts that lie in the interval."
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.max_memory*1000)
    coresMin: $(inputs.cpu)
  - class: DockerRequirement
    dockerPull: 'broadinstitute/gatk:4.2.4.1'
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      $( inputs.input_aligned_reads ? 'gatk' : 'echo gatk' )
  - position: 1
    shellQuote: false
    valueFrom: >-
      CollectReadCounts
  - position: 2
    shellQuote: false
    prefix: "-O"
    valueFrom: >-
      ${var pre = inputs.output_prefix ? inputs.output_prefix : inputs.input_aligned_reads ? inputs.input_aligned_reads.nameroot : 'output'; var ext = inputs.output_format ? inputs.output_format.toLowerCase() : 'hdf5'; return pre+'.'+ext}
inputs:
  reference:
    type: 'File?'
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
  input_aligned_reads:
    type: 'File?'
    secondaryFiles: [{ pattern: ".bai", required: false },{ pattern: "^.bai", required: false },{ pattern: ".crai", required: false },{ pattern: "^.crai", required: false }]
    doc: "BAM/SAM/CRAM file containing reads"
    inputBinding:
      position: 2
      prefix: "-I"
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
  output_format:
    type:
      - 'null'
      - type: enum
        name: output_format
        symbols: ["HDF5","TSV"]
    doc: "Output file format."
    inputBinding:
      position: 2
      prefix: "--format"
  output_prefix:
    type: 'string?'
    doc: "String to use as the prefix for the outputs."
  validation_stringency:
    type:
      - 'null'
      - type: enum
        name: validation_stringency
        symbols: ["SILENT","STRICT","LENIENT"]
    doc: "Validation stringency for all SAM/BAM/CRAM/SRA files read by this program. The default stringency value SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded."
    inputBinding:
      position: 2
      prefix: "--read-validation-stringency"
  max_memory:
    type: 'int?'
    default: 14
    doc: "Maximum GB of RAM to allocate for this tool."
  cpu:
    type: 'int?'
    default: 4
    doc: "Number of CPUs to allocate to this task."
outputs:
  output: { type: 'File?', outputBinding: { glob: "*.hdf5" } }
