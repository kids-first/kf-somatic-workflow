cwlVersion: v1.0
class: CommandLineTool
id: samtools_faidx
doc: |-
  This tool takes an input fasta and optionally a input index for the input fasta.
  If the index is not provided this tool will generate one.
  Finally the tool will return the input reference file with the index (generated or provided) as a secondaryFile.
requirements:
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9'
  - class: InitialWorkDirRequirement
    listing: [$(inputs.input_fasta),$(inputs.input_index)]
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
  - class: ShellCommandRequirement
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      $(inputs.input_index ? 'echo samtools faidx' : 'samtools faidx' )
inputs:
  input_fasta: { type: File, inputBinding: { position: 1 }, doc: "Input fasta file" }
  input_index: { type: 'File?', doc: "Input fasta index" }
outputs:
  fai: 
    type: File
    outputBinding:
      glob: "*.fai" 
