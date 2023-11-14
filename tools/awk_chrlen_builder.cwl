cwlVersion: v1.2
class: CommandLineTool
id: awk_chrlen_builder 
doc: >-
  This tool is an awk command that does the following:
  - Create a list of chromosome names from the input_intervals
  - Print any matching chromosome names from the reference_fai
requirements:
  - class: DockerRequirement
    dockerPull: 'ubuntu:22.04'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
    ramMin: $(inputs.ram * 1000)
  - class: ShellCommandRequirement

baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: |
      awk -F'\t' 'NR==FNR && $1 !~ /^@/ { out[$1]=1; next } { if (out[$1]==1) print $0 }' $(inputs.input_intervals.path) $(inputs.reference_fai.path) > chrlen.tsv

inputs:
  input_intervals: { type: 'File', doc: "Long Reads BAM/CRAM/SAM file."}
  reference_fai: { type: 'File', doc: "FAI index for the FASTA on which the intervals are based." }
  cpu: { type: 'int?', default: 8 }
  ram: { type: 'int?', default: 16 }

outputs:
  chrlen:
    type: File
    outputBinding:
      glob: "chrlen.tsv"
