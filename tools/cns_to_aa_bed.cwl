cwlVersion: v1.2
class: CommandLineTool
id: cns-to-aa-bed
doc: "Basic tool to convert the raw .cns output from CNVkit for use in prepare_aa"
requirements: 
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement 
  - class: DockerRequirement
    dockerPull: 'jluebeck/prepareaa:v0.1203.10'
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4
baseCommand: [/home/programs/PrepareAA-master/scripts/convert_cns_to_bed.py]
inputs:
  input_cns: { type: File, doc: "CNVkit raw cns file to convert", inputBinding: { position: 1, prefix: "--cns_file"} }
outputs:
  aa_cns_bed: 
    type: File
    outputBinding:
      glob: '*_ESTIMATED_PLOIDY_CORRECTED_CN.bed'