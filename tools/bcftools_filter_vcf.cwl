cwlVersion: v1.0
class: CommandLineTool
id: bcftools_filter_vcf
doc: "More generic tool to take in an include expression and optionally an exclude expresssion to filter a vcf"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/bvcftools:latest'
  - class: ResourceRequirement
    ramMin: 1000
    coresMin: 1
  - class: InlineJavascriptRequirement
baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      ${
        var cmd = "bcftools view ";
        if (inputs.include_expression != null){
            cmd += "--include '" + inputs.include_expression + "' " + inputs.input_vcf.path;
            if (inputs.exclude_expression != null){
                cmd += " | bcftools view --exclude '" + inputs.exclude_expression + "' -O z > " + inputs.output_basename + ".vcf.gz;";
            } else {
                cmd += " -O z > " + inputs.output_basename + ".vcf.gz;";
            }
        } else if (inputs.include_expression == null && inputs.exclude_expression != null){
            cmd += "--exclude '" + inputs.exclude_expression + "' " + inputs.input_vcf.path + " -O z > " + inputs.output_basename + ".vcf.gz;";
        } else if (inputs.include_expression == null && inputs.exclude_expression == null){
            cmd = "cp " + inputs.input_vcf.path + " ./" + inputs.output_basename + ".vcf.gz;";
        }
        return cmd;
      }

      tabix $(inputs.output_basename).vcf.gz

inputs:
  input_vcf: File
  include_expression: ['null', string]
  exclude_expression: ['null', string]
  output_basename: string
outputs:
  filtered_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: [.tbi]
