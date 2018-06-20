cwlVersion: v1.0
class: CommandLineTool
id: lancet
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 20000
    coresMin: 16
  - class: DockerRequirement
    dockerPull: 'seandavi/lancet'

baseCommand: [lancet]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      --tumor $(inputs.input_tumor_bam.path)
      --normal $(inputs.input_normal_bam.path)
      --ref $(inputs.reference.path)
      --bed $(inputs.bed.path)
      --num-threads 16 > $(inputs.output_basename).vcf

inputs:
    reference: {type: File, secondaryFiles: [^.dict, .fai]}
    input_tumor_bam: {type: File, secondaryFiles: [.bai]}
    input_normal_bam: {type: File, secondaryFiles: [.bai]}
    bed: {type: File}
    output_basename: {type: string}
outputs:
  - id: output
    type: File
    outputBinding:
      glob: '*.vcf'
