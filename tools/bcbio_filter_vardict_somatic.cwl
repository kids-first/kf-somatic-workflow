cwlVersion: v1.0
class: CommandLineTool
id: bcbio_vardict_fp_somatic_filter
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4
  - class: DockerRequirement
    dockerPull: 'kfdrc/bcbio_vardict_filter'

baseCommand: [python]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      /bcbio_vardict_filter.py $(inputs.input_vcf.path)
      | grep -E "^#|STATUS=StrongSomatic"
      > $(inputs.output_basename).bcbio_vardict_fp_somatic_filter.vcf

      bgzip $(inputs.output_basename).bcbio_vardict_fp_somatic_filter.vcf
      && tabix $(inputs.output_basename).bcbio_vardict_fp_somatic_filter.vcf.gz
inputs:
    input_vcf: File
    output_basename: string

outputs:
  filtered_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: ['.tbi']
