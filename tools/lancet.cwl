cwlVersion: v1.0
class: CommandLineTool
id: lancet
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: ${return inputs.ram * 1000}
    coresMin: 6
  - class: DockerRequirement
    dockerPull: 'kfdrc/lancet:1.0.7'

baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail

      /lancet-1.0.7/lancet
      --tumor $(inputs.input_tumor_bam.path)
      --normal $(inputs.input_normal_bam.path)
      --ref $(inputs.reference.path)
      --bed $(inputs.bed.path)
      --num-threads 6
      --window-size $(inputs.window)
      --padding $(inputs.padding)
      --max-indel-len 50
      > $(inputs.input_tumor_bam.nameroot).$(inputs.bed.nameroot).vcf
      || (echo 'active region filter failed, trying without' && /lancet-1.0.7/lancet
      --tumor $(inputs.input_tumor_bam.path)
      --normal $(inputs.input_normal_bam.path)
      --ref $(inputs.reference.path)
      --bed $(inputs.bed.path)
      --num-threads 6
      --window-size $(inputs.window)
      --active-region-off
      --padding $(inputs.padding)
      --max-indel-len 50
      > $(inputs.input_tumor_bam.nameroot).$(inputs.bed.nameroot).vcf)

inputs:
    reference: {type: File, secondaryFiles: [^.dict, .fai]}
    input_tumor_bam: {type: File, secondaryFiles: [^.bai]}
    input_normal_bam: {type: File, secondaryFiles: [^.bai]}
    ram: {type: ['null', int], default: 12, doc: "Adjust in rare circumstances in which 12 GB is not enough"}
    bed: {type: File}
    output_basename: {type: string}
    window: {type: int, doc: "window size for lancet.  default is 600, recommend 500 for WGS, 600 for exome+"}
    padding: {type: int, doc: "If WGS (less likely), recommend 25, if exome+, recommend half window size"}
outputs:
  lancet_vcf:
    type: File
    outputBinding:
      glob: '*.vcf'
