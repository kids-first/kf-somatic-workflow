cwlVersion: v1.0
class: CommandLineTool
id: samtools_calmd
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/samtools:1.9'
  - class: ResourceRequirement
    ramMin: 10000
    coresMin: 36
  - class: InlineJavascriptRequirement
baseCommand: [samtools, calmd]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -@ $(inputs.threads)
      -b $(inputs.input_reads.path)
      --reference $(inputs.reference.path)
      > $(inputs.input_reads.nameroot).calmd.bam
      && samtools index $(inputs.input_reads.nameroot).calmd.bam $(inputs.input_reads.nameroot).calmd.bai
inputs:
  input_reads: File
  threads: int
  reference: File
outputs:
  bam_file:
    type: File
    outputBinding:
      glob: '*.bam'
    secondaryFiles: [^.bai]
