cwlVersion: v1.2
class: CommandLineTool
id: bedops_catch_megarow
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/bedops:2.4.36'
  - class: ResourceRequirement
    ramMin: 16000
    coresMin: 8
  - class: InlineJavascriptRequirement
baseCommand: [/bin/bash, -c]
arguments:
  - position: 1
    shellQuote: true
    valueFrom: >-
      bgzip -dc $(inputs.in_vcf.path) > temp.vcf
      && vcf2bed --insertions < temp.vcf > temp.ins.vcf
      && vcf2bed --deletions < temp.vcf > temp.del.vcf
      && vcf2bed --snvs < temp.vcf > temp.snv.vcf
inputs:
  in_vcf: {type: File, doc: "VCF to analyse for failure"}
outputs:
  out_vcfs:
    type: 'File[]?'
    outputBinding:
      glob: 'temp*vcf'
