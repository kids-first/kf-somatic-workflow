cwlVersion: v1.0
class: CommandLineTool
id: gatk4_learn_oritentation_bias
label: GATK Learn Bias
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: ${ return inputs.max_memory * 1000 }
    coresMin: 2
baseCommand: [/gatk, LearnReadOrientationModel]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m"
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
  max_memory: {type: 'int?', default: 4, doc: "Maximum memory in GB for GATK LearnReadOrientationModel"}
outputs:
  f1r2_bias:
    type: File
    outputBinding:
      glob: '*.f1r2_bias.tar.gz'
