cwlVersion: v1.0
class: CommandLineTool
id: add_basename 
requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 1000
    coresMin: 1
  - class: ShellCommandRequirement
baseCommand: []

arguments:
  - position: 1
    shellQuote: false
    valueFrom: |
      mv $(inputs.input_file.path) $(inputs.output_basename).$(inputs.input_file.basename)

inputs:
  input_file: { type: File, doc: "Input File" }
  output_basename: { type: string }
outputs:
  output: 
    type: File
    outputBinding:
      glob: $('*' + inputs.input_file.basename)
