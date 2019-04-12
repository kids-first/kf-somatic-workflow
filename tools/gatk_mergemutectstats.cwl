cwlVersion: v1.0
class: CommandLineTool
id: gatk4_mergepileup
label: GATK Merge Stats
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: 4000
    coresMin: 2
baseCommand: [/gatk, MergeMutectStats]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      --java-options "-Xmx3000m"
      -O $(inputs.output_basename).Mutect2.merged.stats 

inputs:
  input_stats:
    type:
      type: array
      items: File
      inputBinding:
        prefix: --stats
    inputBinding:
      position: 1
  output_basename: string
outputs:
  merged_stats:
    type: File
    outputBinding:
      glob: '*.merged.stats'
