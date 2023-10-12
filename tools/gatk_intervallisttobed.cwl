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
baseCommand: [gatk, IntervalListToBed]
inputs:
  input_intervallist: { type: 'File', inputBinding: { position: 2, prefix: "--INPUT" }, doc: "Input IntervalList file." }
  output_filename: { type: 'string?', default: "out.bed", inputBinding: { position: 2, prefix: "--OUTPUT" }, doc: "The output Picard Interval List." }
  score: { type: 'int?', inputBinding: { position: 2, prefix: "--SCORE" }, doc: "The score, between 0-1000, to output for each interval in the BED file." }
  sort: { type: 'boolean?', default: true, inputBinding: { position: 2, prefix: "--SORT", valueFrom: '$(self ? "true" : "false")' }, doc: "If true, sort the output interval list before writing it." }
  ram: { type: 'int?', default: 4, doc: "GB of RAM to allocate to the task." }
  cpu: { type: 'int?', default: 2, doc: "Minimum reserved number of CPU cores for the task." }
outputs:
  output: { type: 'File', outputBinding: { glob: $(inputs.output_filename) } }
