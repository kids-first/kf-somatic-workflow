cwlVersion: v1.0
class: CommandLineTool
id: samtools_cram2bam_plus_calmd
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/samtools:1.9'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 12000
    coresMin: $(inputs.threads)
  
baseCommand: [samtools, view]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -@ $(inputs.threads)
      -h
      -T $(inputs.reference.path)
      $(inputs.input_reads.path)
      | samtools calmd
      -@ 16
      -b
      --reference $(inputs.reference.path)
      -
      > $(inputs.input_reads.nameroot).bam
      && samtools index $(inputs.input_reads.nameroot).bam $(inputs.input_reads.nameroot).bai
inputs:
  input_reads: File
  threads:
    type: ['null', int]
    default: 16
  reference: {type: File, secondaryFiles: [.fai]}
outputs:
  bam_file:
    type: File
    outputBinding:
      glob: '*.bam'
    secondaryFiles: [^.bai]
