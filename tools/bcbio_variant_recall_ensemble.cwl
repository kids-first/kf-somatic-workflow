cwlVersion: v1.0
class: CommandLineTool
id: bcbio_vr
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 6000
    coresMin: 6
  - class: DockerRequirement
    dockerPull: 'migbro/bcbio_vr:0.2.4'

baseCommand: [/bcbio-variation-recall]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      ensemble -c 8
      --names $(inputs.tool_name_csv) $(inputs.output_basename).caller_consensus.vcf.gz
      $(inputs.reference.path)
inputs:
    reference: File
    tool_name_csv: string
    output_basename: string
outputs:
  consensus_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: ['.tbi']
