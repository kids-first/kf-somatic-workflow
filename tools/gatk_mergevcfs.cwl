cwlVersion: v1.0
class: CommandLineTool
id: gatk4_mergevcfs
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/gatk:4.0.12.0'
  - class: ResourceRequirement
    ramMin: 2000
baseCommand: [/gatk, MergeVcfs]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      --java-options "-Xmx2000m"
      --TMP_DIR=./TMP
      --CREATE_INDEX=true
      --SEQUENCE_DICTIONARY=$(inputs.reference_dict.path)
      --OUTPUT=$(inputs.output_basename).$(inputs.tool_name).vcf.gz

inputs:
  input_vcfs:
    type:
      type: array
      items: File
      inputBinding:
        prefix: -I
    secondaryFiles: [.tbi]
    inputBinding:
      position: 1
  reference_dict: File
  tool_name: string
  output_basename: string
outputs:
  merged_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: [.tbi]
