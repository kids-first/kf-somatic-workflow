cwlVersion: v1.0
class: CommandLineTool
id: control-freeC-8.7
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 20000
    coresMin: 4
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/bobo823/control-freec-hg38:v8.7'
baseCommand: [tar, -xzvf]
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
      && /FreeC/freec
      -conf file.cfg
      && mv $(inputs.tumor_bam.basename)_ratio.txt $(inputs.output_basename).ratio.txt
      && mv $(inputs.tumor_bam.basename)_CNVs $(inputs.output_basename).CNVs   
inputs:
  tumor_bam: { type: File, secondaryFiles: [.bai] }
  normal_bam: { type: File, secondaryFiles: [.bai] }
  ref_chrs: File
  chr_len: File
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