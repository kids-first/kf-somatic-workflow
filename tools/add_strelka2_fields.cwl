cwlVersion: v1.0
class: CommandLineTool
id: add_strelka2_fields
doc: >-
  This tool will run a Python script that takes a Strelka2 tumor-normal VCF and adds
      canonical fields like GT and AD
requirements:
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/add-strelka2-fields:1.0.0'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: ${ return inputs.ram * 1000 }
    coresMin: $(inputs.cores)
  - class: ShellCommandRequirement

baseCommand: [/usr/bin/add_strelka2_fields.py]

inputs:
  strelka2_vcf:
    type: File
    inputBinding:
      position: 1
      prefix: '--strelka2_vcf'
    secondaryFiles: .tbi
  tumor_name:
    type: string
    inputBinding:
      position: 2
      prefix: '--tumor_name'
  normal_name:
    type: string
    inputBinding:
      position: 3
      prefix: '--normal_name'
  output_basename:
    type: string
    inputBinding:
      position: 4
      prefix: '--output_basename'
  cores:
    type: int?
    default: 4
  ram:
    type: int?
    default: 3 
    doc: 'RAM requirement in GB'

outputs:
  output:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: .tbi
