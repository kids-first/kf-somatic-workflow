cwlVersion: v1.0
class: CommandLineTool
id: gatk_variantfilter_vcf
doc: "Simple tool add a FILTER tag based on criteria"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0'

baseCommand: [/gatk, VariantFiltration]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      -R $(inputs.reference.path)
      -V $(inputs.input_vcf.path)
      -O $(inputs.output_basename).$(inputs.tool_name).gatk.soft_filtered.vcf.gz
      --filter-name "NORM_HIGH" --filter-expression 'vc.getGenotype("BS_9H6Z0MEG").getAD().1 >0' --filter-name GNOMAD_AF_HIGH --filter-expression "AF > 0.01"
inputs:
    input_vcf: {type: File, secondaryFiles: ['.tbi']}
    reference: {type: File, secondaryFiles: [^.dict]}
    filter_name: 
      type:
        type: array
        items: string
        inputBinding:
          prefix: --filter-name
      inputBinding:
        position: 1
      doc: "Array of names for each filter tag to add"
    filter_expression:
      type:
        type: array
        items: string
        inputBinding:
          prefix: --filter-expression
      inputBinding:
        position: 1
      doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration for clues"
    threads: {type: int, doc: "Number of compression/decompression threads", default: 4}
    output_basename: string
    tool_name: string

outputs:
  gatk_soft_filtered_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: ['.tbi']
