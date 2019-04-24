cwlVersion: v1.0
class: CommandLineTool
id: lancet
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 12000
    coresMin: 6
  - class: DockerRequirement
    dockerPull: 'kfdrc/lancet:1.0.7'

baseCommand: [/lancet-1.0.7/lancet]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      --tumor $(inputs.input_tumor_bam.path)
      --normal $(inputs.input_normal_bam.path)
      --ref $(inputs.reference.path)
      --bed $(inputs.bed.path)
      --active-region-off
      -w 300
      --num-threads 6 >  $(inputs.input_tumor_bam.nameroot).$(inputs.bed.nameroot).vcf

inputs:
    reference: {type: File, secondaryFiles: [^.dict, .fai]}
    input_tumor_bam: {type: File, secondaryFiles: [^.bai]}
    input_normal_bam: {type: File, secondaryFiles: [^.bai]}
    bed: {type: File}
    output_basename: {type: string}
outputs:
  lancet_vcf:
    type: File
    outputBinding:
      glob: '*.vcf'
