cwlVersion: v1.0
class: CommandLineTool
id: gatk4_mergevcfs
label: GATK Merge VCF
doc: "Merge input vcfs"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: 6000
    coresMin: 4
baseCommand: [/gatk, SortVcf]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      --java-options "-Xmx6g"
      -O $(inputs.output_basename).$(inputs.tool_name).merged.vcf
      --SEQUENCE_DICTIONARY $(inputs.reference_dict.path)
      --CREATE_INDEX false

  - position: 2
    shellQuote: false
    valueFrom: >-
      
      && cat $(inputs.output_basename).$(inputs.tool_name).merged.vcf | uniq
      | bgzip > $(inputs.output_basename).$(inputs.tool_name).merged.vcf.gz
      && tabix $(inputs.output_basename).$(inputs.tool_name).merged.vcf.gz

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
      glob: '*.merged.vcf.gz'
    secondaryFiles: [.tbi]
