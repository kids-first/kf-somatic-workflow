cwlVersion: v1.0
class: CommandLineTool
id: gatk4_selectvariants
label: GATK Select PASS
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: 4000
    coresMin: 2
baseCommand: [/gatk, SelectVariants]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      --java-options "-Xmx4000m"
      -V $(inputs.input_vcf.path)
      -O $(inputs.output_basename).$(inputs.tool_name).PASS.vcf.gz
      --select 'vc.isNotFiltered()'

inputs:
  input_vcf: {type: File, secondaryFiles: [.tbi]}
  
outputs:  
  pass_vcf:
    type: File
    outputBinding:
      glob: '*.PASS.vcf.gz'
    secondaryFiles: ['.tbi']

