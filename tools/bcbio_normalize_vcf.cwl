cwlVersion: v1.0
class: CommandLineTool
id: bcbio_normalize_vcf
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/vcfutils:latest'

baseCommand: [bcftools, norm, -m]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      '-any'
      $(inputs.input_vcf.path)
      > $(inputs.output_basename).$(inputs.tool_name).bcf_norm.vcf
      && /vcflib/bin/vcfallelicprimitives -t DECOMPOSED --keep-geno $(inputs.output_basename).$(inputs.tool_name).bcf_norm.vcf
      > $(inputs.output_basename).$(inputs.tool_name).bcf_norm.vcflib_decomp.vcf
      && /vt/vt normalize -r $(inputs.indexed_reference_fasta.path) $(inputs.output_basename).$(inputs.tool_name).bcf_norm.vcflib_decomp.vcf > $(inputs.output_basename).$(inputs.tool_name).bcbio_norm.vcf
      && bgzip $(inputs.output_basename).$(inputs.tool_name).bcbio_norm.vcf
      && tabix $(inputs.output_basename).$(inputs.tool_name).bcbio_norm.vcf.gz
inputs:
    input_vcf: {type: File, secondaryFiles: ['.tbi']}
    indexed_reference_fasta: {type: File, secondaryFiles: ['.fai']}
    output_basename: string
    tool_name: string

outputs:
  bcbio_normalize_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: ['.tbi']
