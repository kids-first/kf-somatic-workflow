cwlVersion: v1.0
class: CommandLineTool
id: control-freeC-R
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 4000
    coresMin: 2
    coresMax: 8
  - class: DockerRequirement
    dockerPull: 'kfdrc/controlfreec:11.5'
baseCommand: [cat]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      /FREEC-11.5/scripts/assess_significance.R
      | R
      --slave
      --args
      $(inputs.cnv_result.path)
      $(inputs.cnv_bam_ratio.path)
      && mv $(inputs.cnv_result.path).p.value.txt $(inputs.cnv_result.basename).p.value.txt
inputs:
  cnv_bam_ratio: File
  cnv_result: File
outputs:
  output_pval:
    type: File
    outputBinding:
      glob: '*.value.txt'