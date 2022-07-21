cwlVersion: v1.2
id: untar_gzip
requirements:
  - class: DockerRequirement
    dockerPull: 'ubuntu:20.04' 
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 4000
    coresMin: 4
class: CommandLineTool
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      mkdir $(inputs.out_dir_name)
      && tar
  - position: 2
    shellQuote: false
    valueFrom: >-
      -xzf
inputs:
  input_tar: { type: File, doc: "Input gzipped tarball to unpack", inputBinding: { position: 2 } }
  out_dir_name: { type: 'string?', doc: "name of output dir", default: "untar_output", inputBinding: { position: 1, prefix: "-C", shellQuote: false} }
  
outputs:
  output:
    type: 'Directory?'
    outputBinding:
      glob: $(inputs.out_dir_name)
