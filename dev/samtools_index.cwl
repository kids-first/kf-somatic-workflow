cwlVersion: v1.0
class: CommandLineTool
id: samtools_index
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/samtools:1.9'
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4
  - class: InlineJavascriptRequirement
baseCommand: [samtools, index]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -@ $(inputs.threads) $(inputs.input_reads.path) $(inputs.input_reads.nameroot).bai
inputs:
  input_reads: File
  threads: int
outputs:
  bai_file:
    type: File
    outputBinding:
      glob: '*.bai'
