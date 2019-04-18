cwlVersion: v1.0
class: CommandLineTool
id: control-freeC-visualize
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 2000
    coresMin: 2
    coresMax: 8
  - class: DockerRequirement
    dockerPull: 'kfdrc/controlfreec:11.5'
baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail

      RATIO=$(inputs.cnv_bam_ratio.path)
      
      if file --mime-type $(inputs.cnv_bam_ratio.path) | grep -q gzip$; then
        gunzip  -c $(inputs.cnv_bam_ratio.path) > ./$(inputs.cnv_bam_ratio.nameroot);
        RATIO=./$(inputs.cnv_bam_ratio.nameroot);
      fi

      cat /FREEC-11.5/scripts/makeGraph.R
      | R
      --slave
      --args 2
      $RATIO
      && mv $RATIO.png $(inputs.output_basename).png
inputs:
  cnv_bam_ratio: File
  output_basename: string
outputs:
  output_png:
    type: File
    outputBinding:
      glob: "$(inputs.output_basename).png"