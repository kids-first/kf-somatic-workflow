cwlVersion: v1.0
class: CommandLineTool
id: tabix_index 
doc: >-
  This tool will run tabix conditionally dependent on whether an index is provided.
  The tool will output the input_file with the index, provided or created within, as a secondary file.
requirements:
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9'
  - class: InitialWorkDirRequirement
    listing: [$(inputs.input_file),$(inputs.input_index)]
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
  - class: ShellCommandRequirement

baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      $(inputs.input_file ? inputs.input_index ? 'echo tabix -p vcf' : 'tabix -p vcf' : 'echo tabix -p vcf')

inputs:
  input_file: { type: 'File?', doc: "Position sorted and compressed by bgzip input file", inputBinding: { position: 1, shellQuote: false } }
  input_index: { type: 'File?', doc: "Index file for the input_file, if one exists" }

outputs:
  output: 
    type: File?
    outputBinding:
      glob: "*.gz"
    secondaryFiles: [.tbi]
