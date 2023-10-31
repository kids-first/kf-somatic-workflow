cwlVersion: v1.2
class: CommandLineTool
id: gatk_bedtointervallist
doc: |
  Converts a BED file to a Picard Interval List.  This tool provides easy
  conversion from BED to the Picard interval_list format which is required by
  many Picard processing tools. Note that the coordinate system of BED files is
  such that the first base or position in a sequence is numbered "0", while in
  interval_list files it is numbered "1"
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.max_memory * 1000)
    coresMin: $(inputs.cores)
  - class: DockerRequirement
    dockerPull: 'broadinstitute/gatk:4.4.0.0'
baseCommand: [gatk, BedToIntervalList]
inputs:
  input_bed: { type: 'File', inputBinding: { position: 2, prefix: "--INPUT" }, doc: "The input BED file." }
  output_filename: { type: 'string?', default: "out.interval_list", inputBinding: { position: 2, prefix: "--OUTPUT" }, doc: "The output Picard Interval List." }
  reference_dict: { type: 'File', inputBinding: { position: 2, prefix: "--SEQUENCE_DICTIONARY" }, doc: "The sequence dictionary, or BAM/VCF/IntervalList from which a dictionary can be extracted." }
  sort: { type: 'boolean?', default: true, inputBinding: { position: 2, prefix: "--SORT", valueFrom: '$(self ? "true" : "false")' }, doc: "If true, sort the output interval list before writing it." }
  unique: { type: 'boolean?', default: false, inputBinding: { position: 2, prefix: "--UNIQUE", valueFrom: '$(self ? "true" : "false")' }, doc: "If true, unique the output interval list by merging overlapping regions, before writing it (implies sort=true)." }
  max_memory: { type: 'int?', default: 4, doc: "GB of RAM to allocate to the task." }
  cores: { type: 'int?', default: 2, doc: "Minimum reserved number of CPU cores for the task." }
outputs:
  output: { type: 'File', outputBinding: { glob: $(inputs.output_filename) } }
