cwlVersion: v1.0
class: CommandLineTool
id: gatk4_calulcate_contamination
label: GATK Calculate Contamination
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: 4000
    coresMin: 2
baseCommand: [/gatk, CalculateContamination]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      --java-options "-Xmx4000m"
      -I $(inputs.tumor_pileup.path)
      --matched-normal $(inputs.normal_pileup.path)
      -O $(inputs.output_basename).contamination.table
      --tumor-segmentation $(inputs.output_basename).segmentation.table

inputs:
  tumor_pileup: File
  normal_pileup: File
  output_basename: string
outputs:
  contamination_table:
    type: File
    outputBinding:
      glob: '*.contamination.table'
  segmentation_table:
    type: File
    outputBinding:
      glob: '*.segmentation.table'

