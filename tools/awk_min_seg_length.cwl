cwlVersion: v1.2
id: awk_min_seg_length
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'ubuntu:20.04'
class: CommandLineTool
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      awk -F'\t' -v min=$(inputs.default_min_len) '$1 ~ /chr/ && $3-$2 < min {min = $3-$2} END {print min}' $(inputs.input_file.path) > value.txt
inputs:
  input_file: { type: 'File' }
  default_min_len: { type: 'int?', default: 150 }
outputs:
  output:
    type: 'int'
    outputBinding:
      glob: 'value.txt'
      loadContents: true
      outputEval: $(parseInt(self[0].contents))
