cwlVersion: v1.0
class: CommandLineTool
id: samtools_cram2bam
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/samtools:1.9'
  - class: ResourceRequirement
    ramMin: 10000
    coresMin: 36
  - class: InlineJavascriptRequirement
baseCommand: [samtools, view]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -@ $(inputs.threads)
      -bh $(inputs.input_reads.path)
      -m 1G
      -T $(inputs.reference.path)
      > $(inputs.input_reads.nameroot).bam
      && samtools index $(inputs.input_reads.nameroot).bam
inputs:
  input_reads: File
  threads: int
  reference: File
outputs:
  bam_file:
    type: File
    outputBinding:
      glob: '*.bam'
    secondaryFiles: [.bai]
