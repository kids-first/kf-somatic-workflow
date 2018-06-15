cwlVersion: v1.0
class: CommandLineTool
id: gatk4_Mutect2
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/gatk:4.0.3.0'
  - class: ResourceRequirement
    ramMin: 8000
baseCommand: [/gatk, Mutect2]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      --java-options "-Xms8000m"
      -R $(inputs.reference.path)
      -I $(inputs.input_tumor_bam.path)
      -I $(inputs.input_normal_bam.path)
      -tumor $(inputs.input_tumor_name)
      -normal $(inputs.input_normal_name)
      --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter
      -L $(inputs.interval_list.path)
      -O $(inputs.input_tumor_bam.nameroot).vcf.gz
      -bamout $(inputs.input_tumor_bam.nameroot)_tumor_normal.bam

inputs:
  reference: {type: File, secondaryFiles: [^.dict, .fai]}
  input_tumor_bam: {type: File, secondaryFiles: [^.bai]}
  input_tumor_name: string
  input_normal_bam: {type: File, secondaryFiles: [^.bai]}
  input_normal_name: string
  interval_list: File

outputs:
  output:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: [.tbi]
