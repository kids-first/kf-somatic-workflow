cwlVersion: v1.0
class: CommandLineTool
id: samtools_calmd
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/samtools:1.9'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 32000
    coresMin: inputs.threads
baseCommand: [samtools, calmd]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -@ 16
      -b $(inputs.input_reads.path)
      --reference $(inputs.reference.path)
      > $(inputs.input_reads.nameroot).calmd.bam
      && samtools index $(inputs.input_reads.nameroot).calmd.bam $(inputs.input_reads.nameroot).calmd.bai
inputs:
  input_reads: File
  threads:
    type: ['null', int]
    default: 16
  reference: File
outputs:
  bam_file:
    type: File
    outputBinding:
      glob: '*.bam'
    secondaryFiles: [^.bai]
