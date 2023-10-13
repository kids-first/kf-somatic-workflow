cwlVersion: v1.2
class: CommandLineTool
id: gatk_bedtointervallist
doc: |
  Converts an Picard IntervalList file to a BED file.
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.cpu)
  - class: DockerRequirement
    dockerPull: 'broadinstitute/gatk:4.4.0.0'
baseCommand: []
arguments:
- position: 0
  shellQuote: false
  valueFrom: >
    gatk IntervalListToBed
- position: 10
  shellQuote: false
  prefix: '&&'
  valueFrom: >
    bgzip --stdout --threads  $(inputs.cpu)
- position: 19
  shellQuote: false
  valueFrom: >
    $(inputs.output_filename) > $(inputs.output_filename).gz
- position: 20
  shellQuote: false
  prefix: '&&'
  valueFrom: >
    tabix --preset bed
- position: 28
  shellQuote: false
  valueFrom: >
    $(inputs.output_filename).gz

inputs:
  input_intervallist: { type: 'File', inputBinding: { position: 2, prefix: "--INPUT" }, doc: "Input IntervalList file." }
  output_filename: { type: 'string?', default: "out.bed", inputBinding: { position: 2, prefix: "--OUTPUT" }, doc: "The output Picard Interval List." }

  # IntervalListToBed options
  score: { type: 'int?', inputBinding: { position: 2, prefix: "--SCORE" }, doc: "The score, between 0-1000, to output for each interval in the BED file." }
  sort: { type: 'boolean?', default: true, inputBinding: { position: 2, prefix: "--SORT", valueFrom: '$(self ? "true" : "false")' }, doc: "If true, sort the output interval list before writing it." }

  # BGZIP Options
  offset: { type: 'int?', inputBinding: { position: 12, prefix: "--offset"}, doc: "decompress at virtual file pointer (0-based uncompressed offset)" }
  size: { type: 'int?', inputBinding: { position: 12, prefix: "--size"}, doc: "decompress INT bytes (uncompressed size)" }

  # Tabix Options
  zero_based: { type: 'boolean?', inputBinding: { position: 22, prefix: "--zero-based"}, doc: "coordinates are zero-based" }
  begin: { type: 'int?', inputBinding: { position: 22, prefix: "--begin"}, doc: "column number for region start [4]" }
  comment: { type: 'string?', inputBinding: { position: 22, prefix: "--comment"}, doc: "skip comment lines starting with CHAR [null]" }
  csi: { type: 'boolean?', inputBinding: { position: 22, prefix: "--csi"}, doc: "generate CSI index for VCF (default is TBI)" }
  end: { type: 'int?', inputBinding: { position: 22, prefix: "--end"}, doc: "column number for region end (if no end, set INT to -b) [5]" }
  min_shift: { type: 'int?', inputBinding: { position: 22, prefix: "--min-shift"}, doc: "set minimal interval size for CSI indices to 2^INT [14]" }
  sequence: { type: 'int?', inputBinding: { position: 22, prefix: "--sequence"}, doc: "column number for sequence names (suppressed by -p) [1]" }
  skip_lines: { type: 'int?', inputBinding: { position: 22, prefix: "--skip-lines"}, doc: "skip first INT lines [0]" }

  # Resource Control
  ram: { type: 'int?', default: 4, doc: "GB of RAM to allocate to the task." }
  cpu: { type: 'int?', default: 2, doc: "Minimum reserved number of CPU cores for the task." }
outputs:
  outbed: { type: 'File', outputBinding: { glob: $(inputs.output_filename) } }
  outbedgz: { type: 'File', secondaryFiles: [{pattern: '.tbi', required: false}, {pattern: '.csi', required: false}], outputBinding: { glob: $(inputs.output_filename).gz } }
