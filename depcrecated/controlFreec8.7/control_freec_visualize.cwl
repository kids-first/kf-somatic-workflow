cwlVersion: v1.0
class: CommandLineTool
id: control-freeC-visualize
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 10000
    coresMin: 2
    coresMax: 8
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/bobo823/rtools:v3.2.3'
baseCommand: [cat]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      /makeGraph.R
      | R
      --slave
      --args 2
      $(inputs.cnv_bam_ratio.path)
      && mv $(inputs.cnv_bam_ratio.path).png $(inputs.output_basename).png
inputs:
  cnv_bam_ratio: File
  output_basename: string
outputs:
  output_png:
    type: File
    outputBinding:
      glob: '*png'