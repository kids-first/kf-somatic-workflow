cwlVersion: v1.0
class: CommandLineTool
id: bcftools_strip_info
doc: "Quick tool to strip info from vcf file before re-annotation"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/vcfutils:latest'

baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      set -eo pipefail

      (bcftools annotate -x $(inputs.strip_info) $(inputs.input_vcf.path) -O z 
      -o $(inputs.output_basename).$(inputs.tool_name).INFO_stripped.vcf.gz &&
      tabix $(inputs.output_basename).$(inputs.tool_name).INFO_stripped.vcf.gz) ||
      (echo "Check errors, likely does not have INFO, trying to pass input instead" >&2; cp $(inputs.input_vcf.path) .;
      cp $(inputs.input_vcf.secondaryFiles[0].path)  .;)

inputs:
    input_vcf: {type: File, secondaryFiles: ['.tbi']}
    output_basename: string
    tool_name: string
    strip_info: {type: ['null', string], doc: "If given, remove previous annotation information based on INFO file, i.e. to strip VEP info, use INFO/ANN"}

outputs:
  stripped_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: ['.tbi']
