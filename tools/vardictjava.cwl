cwlVersion: v1.0
class: CommandLineTool
id: vardictjava
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 20
    coresMin: 4
  - class: DockerRequirement
    dockerPull: 'migbro/vardict:1.5.8'

baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail

      export VAR_DICT_OPTS='"-Xms768m" "-Xmx18g"' && /VarDict-1.5.8/bin/VarDict
      -G $(inputs.reference.path)
      -f 0.01 -th 8 --nosv
      -N $(inputs.output_basename)
      -b '$(inputs.input_tumor_bam.path)|$(inputs.input_normal_bam.path)'
      -z -c 1 -S 2 -E 3 -g 4 $(inputs.bed.path) -y
      | /VarDict-1.5.8/bin/testsomatic.R
      | /VarDict-1.5.8/bin/var2vcf_paired.pl
      -N '$(inputs.input_tumor_name)|$(inputs.input_normal_name)'
      -f 0.01 > $(inputs.output_basename).vcf &&
      bgzip  $(inputs.output_basename).vcf &&
      tabix  $(inputs.output_basename).vcf.gz

inputs:
  reference: {type: File, secondaryFiles: [^.dict, .fai]}
  input_tumor_bam: {type: File, secondaryFiles: [^.bai]}
  input_tumor_name: string
  input_normal_bam: {type: File, secondaryFiles: [^.bai]}
  input_normal_name: string
  output_basename: {type: string}
  bed: {type: File}
outputs:
  vardict_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: [.tbi]