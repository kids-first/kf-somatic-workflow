cwlVersion: v1.0
class: CommandLineTool
id: bcftools_reheader_vcf
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/bvcftools:latest'
  - class: ResourceRequirement
    ramMin: 1000
    coresMin: 1
  - class: InlineJavascriptRequirement
baseCommand: [echo]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.input_normal_name) > sample_list.txt &&
      echo $(inputs.input_tumor_name) >> sample_list.txt &&
      bcftools reheader -s sample_list.txt $(inputs.input_vcf.path) > $(inputs.input_vcf.nameroot.replace(".vcf", ".reheadered.vcf.gz")) &&
      tabix $(inputs.input_vcf.nameroot.replace(".vcf", ".reheadered.vcf.gz"))

inputs:
  input_vcf: File
  input_normal_name: string
  input_tumor_name: string
  old_tumor_name: { type: 'string?', doc: "Not needed for this tool. Just here to use as a workflow flag"}
outputs:
  reheadered_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: [.tbi]
