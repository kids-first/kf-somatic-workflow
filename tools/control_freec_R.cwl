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

      cat /FREEC-11.5/scripts/assess_significance.R
      | R
      --slave
      --args
      $(inputs.cnv_result.path)
      $RATIO
      && mv $(inputs.cnv_result.path).p.value.txt $(inputs.cnv_result.basename).p.value.txt
inputs:
  cnv_bam_ratio: File
  cnv_result: File
outputs:
  output_pval:
    type: File
    outputBinding:
      glob: '*.value.txt'