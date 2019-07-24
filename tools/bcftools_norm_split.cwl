cwlVersion: v1.0
class: CommandLineTool
id: bcftools_norm_split_vcf
doc: "Normalized a vcf by splitting multi-allelic locations, then splits into snps+mnps and indels files"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/bvcftools:latest'
  - class: ResourceRequirement
    ramMin: 1000
    coresMin: 1
  - class: InlineJavascriptRequirement
baseCommand: [bcftools,norm,-c,w,-m,-any,-Oz,-f]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.reference.path) $(inputs.input_vcf.path) -Oz >  $(inputs.input_vcf.nameroot.replace('.vcf','')).bcfNorm.vcf.gz &&
      bcftools view --types snps,mnps $(inputs.input_vcf.nameroot.replace('.vcf','')).bcfNorm.vcf.gz -Oz > $(inputs.input_vcf.nameroot.replace('.vcf','')).bcfNorm.snv_mnv.vcf.gz &&
      bcftools view --types indels $(inputs.input_vcf.nameroot.replace('.vcf','')).bcfNorm.vcf.gz -Oz > $(inputs.input_vcf.nameroot.replace('.vcf','')).bcfNorm.indels.vcf.gz &&
      ls *.vcf.gz | xargs -IFN tabix FN


inputs:
  input_vcf: File
  reference: File
outputs:
  normalized_vcf:
    type: File
    outputBinding:
      glob: '*.bcfNorm.vcf.gz'
    secondaryFiles: [.tbi]
  snv_mnv_vcf:
    type: File
    outputBinding:
      glob: 'bcfNorm.snv_mnv.vcf.gz'
  indelvcf:
    type: File
    outputBinding:
      glob: 'bcfNorm.indels.vcf.gz'
