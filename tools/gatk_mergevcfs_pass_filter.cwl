cwlVersion: v1.0
class: CommandLineTool
id: gatk4_mergevcfs_selectvcfs
label: GATK Merge and PASS
doc: "Merge input vcfs and output only PASS"
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
      --OUTPUT=$(inputs.output_basename).$(inputs.tool_name).vcf.gz &&
      /gatk SelectVariants
      --java-options "-Xmx2000m"
      -V $(inputs.output_basename).$(inputs.tool_name).vcf.gz
      -O $(inputs.output_basename).$(inputs.tool_name).PASS.vcf.gz
      --select 'vc.isNotFiltered()'

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
      glob: '*PASS.vcf.gz'
    secondaryFiles: [.tbi]