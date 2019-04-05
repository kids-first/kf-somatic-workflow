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
    dockerPull: 'migbro/controlfreec:11.5'

baseCommand: [tar, -xzf]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.ref_chrs.path)
      && perl /get_config.pl
      $(inputs.chr_len.path)
      ./GRCh38_everyChrs
      $(inputs.normal_bam.path)
      $(inputs.tumor_bam.path)
      $(inputs.threads)
      && /FREEC-11.5/src/freec
      -conf file.cfg
      && mv $(inputs.tumor_bam.basename)_ratio.txt $(inputs.output_basename).ratio.txt
      && mv $(inputs.tumor_bam.basename)_CNVs $(inputs.output_basename).CNVs
inputs:
  tumor_bam: { type: File, secondaryFiles: [^.bai] }
  normal_bam: { type: File, secondaryFiles: [^.bai] }
  ref_chrs: File
  chr_len: File
  threads: int
  output_basename: string
outputs:
  output_txt:
    type: File
    outputBinding:
      glob: '*.ratio.txt'
  output_cnv:
    type: File
    outputBinding:
      glob: '*.CNVs'