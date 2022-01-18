cwlVersion: v1.1
class: CommandLineTool
id: gatk_collectreadcounts
doc: >-
  Collects reference and alternate allele counts at specified sites. The alt count is defined as the
  total count minus the ref count, and the alt nucleotide is defined as the non-ref base with the
  highest count, with ties broken by the order of the bases in AllelicCountCollector#BASES. Only
  reads that pass the specified read filters and bases that exceed the specified
  minimum-base-quality will be counted.
  Notes: The input interval list to this file should be either a dbSNP VCF+TBI or the germline VCF+TBI from of the normal inputBAM.
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.max_memory*1000)
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
    prefix: "--java-options"
    valueFrom: >-
      $("\"-Xmx"+Math.floor(inputs.max_memory*1000/1.074 - 1)+"M\"")
  - position: 2
    shellQuote: false
    valueFrom: CollectAllelicCounts
  - position: 3
    shellQuote: false
    prefix: "-O"
    valueFrom: >-
      ${var pre = inputs.output_prefix ? inputs.output_prefix : inputs.input_aligned_reads ? inputs.input_aligned_reads.nameroot : 'output'; return pre+'.allelicCounts.tsv'}
inputs:
  reference:
    type: 'File?'
    doc: "Reference fasta"
    secondaryFiles: ['.fai']
    inputBinding:
      position: 3
      prefix: "-R"
  sequence_dictionary:
    type: 'File?'
    doc: "Use the given sequence dictionary as the master/canonical sequence dictionary. Must be a .dict file."
    inputBinding:
      position: 3
      prefix: "--sequence-dictionary"
  input_aligned_reads:
    type: 'File?'
    secondaryFiles: [{ pattern: ".bai", required: false },{ pattern: "^.bai", required: false },{ pattern: ".crai", required: false },{ pattern: "^.crai", required: false }]
    doc: "BAM/SAM/CRAM file containing reads"
    inputBinding:
      position: 3
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
      position: 3
  input_interval_list:
    type: 'File?'
    secondaryFiles: [{ pattern: ".tbi", required: false }]
    doc: "One or more genomic intervals over which to operate. Use this input when providing interval list files or other file based inputs."
    inputBinding:
      position: 3
      prefix: "-L"
  interval_padding:
    type: 'int?'
    doc: "Amount of padding (in bp) to add to each interval you are including."
    inputBinding:
      position: 3
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
      position: 3
  input_exclude_interval_list:
    type: 'File?'
    secondaryFiles: [{ pattern: ".tbi", required: false }]
    doc: "One or more genomic intervals to exclude from processing. Use this input when providing interval list files or other file based inputs."
    inputBinding:
      position: 3
      prefix: "-XL"
  interval_exclusion_padding:
    type: 'int?'
    doc: "Amount of padding (in bp) to add to each interval you are excluding."
    inputBinding:
      position: 3
      prefix: "-ixp"
  interval_merging_rule:
    type:
      - 'null'
      - type: enum
        name: interval_merging_rule
        symbols: ["ALL","OVERLAPPING_ONLY"]
    doc: "By default, the program merges abutting intervals (i.e. intervals that are directly side-by-side but do not actually overlap) into a single continuous interval. However you can change this behavior if you want them to be treated as separate intervals instead."
    inputBinding:
      position: 3
      prefix: "-imr"
  interval_set_rule:
    type:
      - 'null'
      - type: enum
        name: interval_set_rule
        symbols: ["UNION","INTERSECTION"]
    doc: "By default, the program will take the UNION of all intervals specified using -L and/or -XL. However, you can change this setting for -L, for example if you want to take the INTERSECTION of the sets instead. E.g. to perform the analysis only on chromosome 1 exomes, you could specify -L exomes.intervals -L 1 --interval-set-rule INTERSECTION. However, it is not possible to modify the merging approach for intervals passed using -XL (they will always be merged using UNION). Note that if you specify both -L and -XL, the -XL interval set will be subtracted from the -L interval set."
    inputBinding:
      position: 3
      prefix: "-isr"
  max_depth_per_sample:
    type: 'int?'
    doc: "Maximum number of reads to retain per sample per locus. Reads above this threshold will be downsampled. Set to 0 to disable."
    inputBinding:
      position: 3
      prefix: "--max-depth-per-sample"
  minimum_base_quality:
    type: 'int?'
    doc: "Minimum base quality. Base calls with lower quality will be filtered out of pileups."
    inputBinding:
      position: 3
      prefix: "--minimum-base-quality"
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
      position: 3
      prefix: "--read-validation-stringency"
  max_memory:
    type: 'int?'
    default: 20
    doc: "Maximum GB of RAM to allocate for this tool."
outputs:
  output: { type: 'File?', outputBinding: { glob: "*.allelicCounts.tsv" } }
