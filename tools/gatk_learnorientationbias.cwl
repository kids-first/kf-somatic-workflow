cwlVersion: v1.0
class: CommandLineTool
id: gatk4_learn_oritentation_bias
label: GATK Learn Bias
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: 4000
    coresMin: 2
baseCommand: [/gatk, LearnReadOrientationModel]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      --java-options "-Xmx4000m"
      -O $(inputs.output_basename).$(inputs.tool_name).f1r2_bias.tar.gz 

inputs:
  input_tgz:
    type:
      type: array
      items: File
      inputBinding:
        prefix: -I
    inputBinding:
      position: 1
  tool_name: string
  output_basename: string
outputs:
  f1r2_bias:
    type: File
    outputBinding:
      glob: '*.f1r2_bias.tar.gz'
