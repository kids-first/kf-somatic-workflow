cwlVersion: v1.0
class: CommandLineTool
id: cnvnator
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 8000
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/cnvnator:latest'
baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      ln -s $(inputs.reference.path) ./ &&
      perl /opt/split_fasta.pl $(inputs.reference.path) ./ &&
      /usr/local/bin/cnvnator -root $(inputs.input_bam.nameroot).root -tree $(inputs.input_bam.path) -unique &&
      /usr/local/bin/cnvnator -root $(inputs.input_bam.nameroot).root -his 100 -d ./ && 
      /usr/local/bin/cnvnator -root $(inputs.input_bam.nameroot).root -stat 100   && 
      /usr/local/bin/cnvnator -root $(inputs.input_bam.nameroot).root -partition 100  && 
      /usr/local/bin/cnvnator -root $(inputs.input_bam.nameroot).root -call 100 > $(inputs.output_basename).cnvnator.vcf 
      
inputs:
  reference: {type: File, secondaryFiles: [^.dict, .fai]}
  input_bam: {type: File, secondaryFiles: [.crai]}
  output_basename: string
outputs:
  output:
    type: File
    outputBinding:
      glob: '*.vcf'
