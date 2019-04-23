cwlVersion: v1.0
class: CommandLineTool
id: gatk4_mergevcfs
label: Split bed
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/bvcftools:latest'
  - class: ResourceRequirement
    ramMin: 6000
    coresMin: 4
baseCommand: [gunzip, -c]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      $(inputs.input_strelka_bed.path) | split -l 1 --additional-suffix .bed -

inputs:
  input_strelka_bed: File
outputs:
  output:
    type: 'File[]'
    outputBinding:
      glob: '*.bed'
