cwlVersion: v1.0
class: CommandLineTool
id: vardictjava
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 10000
  - class: DockerRequirement
    dockerPull: 'jinh2/vardictjava'

baseCommand: [/opt/VarDictJava/build/install/VarDict/bin/VarDict]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -G $(inputs.reference.path)
      -f 0.01 -th 8
      -N $(inputs.output_basename)
      -b '$(inputs.input_tumor_bam.path)|$(inputs.input_normal_bam.path)'
      -z -c 1 -S 2 -E 3 -g 4 $(inputs.bed.path)
      | /opt/VarDictJava/VarDict/testsomatic.R
      | /opt/VarDictJava/VarDict/var2vcf_paired.pl
      -N '$(inputs.input_tumor_BSID)|$(inputs.input_normal_BSID)'
      -f 0.01 > $(inputs.output_basename).vcf

inputs:
    reference: {type: File, secondaryFiles: [^.dict, .fai]}
    input_tumor_BSID: {type: string}
    input_normal_BSID: {type: string}
    input_tumor_bam: {type: File, secondaryFiles: [.bai]}
    input_normal_bam: {type: File, secondaryFiles: [.bai]}
    output_basename: {type: string}
    bed: {type: File}
outputs:
  - id: output
    type: File
    outputBinding:
      glob: '*.vcf'
      