cwlVersion: v1.0
class: CommandLineTool
id: control-freeC-11-5

requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 20000
    coresMin: 16
  - class: DockerRequirement
    dockerPull: 'kfdrc/controlfreec:11.5'

baseCommand: [tar, -xzf]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.ref_chrs.path)
      && /FREEC-11.5/src/freec
      -conf $(inputs.config_file.path)
      && mv $(inputs.tumor_bam.basename)_ratio.txt $(inputs.output_basename).ratio.txt
      && gzip $(inputs.output_basename).ratio.txt
      && mv $(inputs.tumor_bam.basename)_CNVs $(inputs.output_basename).CNVs

inputs:
  tumor_bam: { type: File, secondaryFiles: [^.bai] }
  normal_bam: { type: File, secondaryFiles: [^.bai] }
  ref_chrs: File
  chr_len: File
  threads: int
  output_basename: string
  config_file: File
  capture_regions: {type: ['null', File], doc: "If not WGS, provide "}

outputs:
  output_txt:
    type: File
    outputBinding:
      glob: '*.ratio.txt.gz'
  output_cnv:
    type: File
    outputBinding:
      glob: '*.CNVs'