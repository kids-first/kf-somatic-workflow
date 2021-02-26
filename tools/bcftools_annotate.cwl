cwlVersion: v1.0
class: CommandLineTool
id: bcftools_annotate_vcf
doc: "Simple tool to annotate a vcf using bcftools and an annotation vcf"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/vcfutils:latest'

baseCommand: [bcftools, annotate]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      --annotations $(inputs.annotation_vcf.path)
      --columns $(inputs.columns)
      -o $(inputs.output_basename).$(inputs.tool_name).bcf_annotated.vcf.gz
      -O z
      --threads $(inputs.threads)
      $(inputs.input_vcf.path)
      && tabix $(inputs.output_basename).$(inputs.tool_name).bcf_annotated.vcf.gz
inputs:
    input_vcf: {type: File, secondaryFiles: ['.tbi']}
    annotation_vcf: {type: File, secondaryFiles: ['.tbi'], doc: "bgzipped annotation vcf file"}
    columns: {type: string, doc: "csv string of columns from annotation to port into the input vcf, i.e INFO/AF"}
    threads: {type: int?, doc: "Number of compression/decompression threads", default: 4}
    output_basename: string
    tool_name: string

outputs:
  bcftools_annotated_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: ['.tbi']
