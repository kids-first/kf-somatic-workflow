class: CommandLineTool
cwlVersion: v1.0
id: bwa_index
doc: |-
  This tool conditionally generates the bwa 64 indexes for an input fasta file using bwa index.
  The tool will generate the indexes only of generate_bwa_indexes is set to true AND any of the alt,
  amb, ann, bwt, pac, or sa files is missing.
  The tool returns the six indexes as its output.
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
  - class: InitialWorkDirRequirement
    listing: [$(inputs.input_fasta),$(inputs.input_alt),$(inputs.input_amb),$(inputs.input_ann),$(inputs.input_bwt),$(inputs.input_pac),$(inputs.input_sa)]
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/bwa:0.7.17-dev'
  - class: InlineJavascriptRequirement
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      $(inputs.input_alt && inputs.input_amb && inputs.input_ann && inputs.input_bwt && inputs.input_pac && inputs.input_sa ? 'echo bwa' : inputs.generate_bwa_indexes ? 'bwa' : 'echo bwa')
  - position: 1
    shellQuote: false
    valueFrom: >-
      index -6 -a bwtsw 
inputs:
  generate_bwa_indexes: { type: 'boolean?' }
  input_fasta: { type: File, inputBinding: { position: 2, valueFrom: $(self.basename) } }
  input_alt: { type: 'File?' }
  input_amb: { type: 'File?' }
  input_ann: { type: 'File?' }
  input_bwt: { type: 'File?' }
  input_pac: { type: 'File?' }
  input_sa: { type: 'File?' }

outputs:
  alt:
    type: 'File?'
    outputBinding:
      glob: '*.64.alt'
  amb:
    type: 'File?'
    outputBinding:
      glob: '*.64.amb'
  ann:
    type: 'File?'
    outputBinding:
      glob: '*.64.ann'
  bwt:
    type: 'File?'
    outputBinding:
      glob: '*.64.bwt'
  pac:
    type: 'File?'
    outputBinding:
      glob: '*.64.pac'
  sa:
    type: 'File?'
    outputBinding:
      glob: '*.64.sa'
