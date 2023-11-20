cwlVersion: v1.0
class: CommandLineTool
id: gatk4_selectvariants
label: GATK Select PASS
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4
  - class: InitialWorkDirRequirement
    listing:
      - entryname: select_pass.sh
        entry: |
          #!/usr/bin/env bash

          set -eo pipefail

          if [[ $(inputs.mode) = gatk ]]
          then
            /gatk SelectVariants --java-options "-Xmx8000m" --exclude-filtered TRUE -V $(inputs.input_vcf.path) -O $(inputs.output_basename).$(inputs.tool_name).PASS.vcf.gz
          elif [[ $(inputs.mode) = grep ]]
          then
            zcat $(inputs.input_vcf.path) | grep -E "^#|PASS" | bgzip > $(inputs.output_basename).$(inputs.tool_name).PASS.vcf.gz && tabix $(inputs.output_basename).$(inputs.tool_name).PASS.vcf.gz
          else
            echo "$(inputs.mode) is not a valid mode. Choices are gatk or grep" >&2 && exit 1
          fi

baseCommand: ["/bin/bash", "-x", "select_pass.sh"]

inputs:
  input_vcf: {type: File, secondaryFiles: [.tbi]}
  output_basename: string
  tool_name: string
  mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression", default: "gatk"}
  
outputs:  
  pass_vcf:
    type: File
    outputBinding:
      glob: '*.PASS.vcf.gz'
    secondaryFiles: ['.tbi']

