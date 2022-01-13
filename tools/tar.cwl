cwlVersion: v1.0
requirements:
  - class: DockerRequirement
    dockerPull: 'ubuntu:18.04' 
  - class: InlineJavascriptRequirement
class: CommandLineTool
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      $( inputs.output_file ? 'tar' : 'echo tar' )
inputs:
  create:
    type: 'boolean?'
    default: true
    inputBinding:
      prefix: --create
      position: 1
  output_file:
    type: 'string?'
    inputBinding:
      prefix: --file
      position: 1
  gzip:
    type: 'boolean?'
    default: true
    inputBinding:
      prefix: --gzip
      position: 1
  input_files:
    type: 'File[]?'
    inputBinding:
      position: 2
outputs:
  output:
    type: 'File?'
    outputBinding:
      glob: $(inputs.output_file)
