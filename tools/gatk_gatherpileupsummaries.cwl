cwlVersion: v1.0
class: CommandLineTool
id: gatk4_mergepileup
label: GATK Merge Pileups
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: 4000
    coresMin: 2
baseCommand: [/gatk, GatherPileupSummaries]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      --java-options "-Xmx3000m"
      --sequence-dictionary $(inputs.reference_dict.path)
      -O $(inputs.output_basename).$(inputs.tool_name).merged.pileup.table 

inputs:
  input_tables:
    type:
      type: array
      items: File
      inputBinding:
        prefix: -I
    inputBinding:
      position: 1
  reference_dict: File
  tool_name: string
  output_basename: string
outputs:
  merged_table:
    type: File
    outputBinding:
      glob: '*.merged.pileup.table'
