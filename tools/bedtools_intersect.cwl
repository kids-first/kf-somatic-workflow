cwlVersion: v1.0
class: CommandLineTool
id: bedtools_intersect
doc: "Intersect VCF with bedtools intersect, i.e. for WXS filter germline on unpadded intervals to ensure accurate CNV calls"
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

      ${
          var cmd = " >&2 echo checking if skip intersect flag given;";
          var out_vcf = inputs.output_basename + ".bed_intersect.vcf";
          if (inputs.input_vcf == null){
            return "echo No vcf provided, skipping >&2 && exit 0;"
          }
          if (inputs.flag != "N"){
            cmd += "bedtools intersect -a " + inputs.input_vcf.path + " -b " + inputs.input_bed_file.path + " -header -wa > " + out_vcf + " && ";
            cmd +=  "bgzip " + out_vcf + " && tabix " + out_vcf + ".gz;";
          }else{
            cmd += " >&2 echo Value N given, passing through vcf; cp " + inputs.input_vcf.path + " " + out_vcf + ".gz;"
            cmd += "tabix " + out_vcf + ".gz;"
          }
          return cmd;
      }

inputs:
    input_vcf: {type: 'File?', secondaryFiles: ['.tbi'], doc: "Input VCF file. Skips if not provided."}
    input_bed_file: File
    output_basename: string
    flag: {type: ['null', string], doc: "If N, skip and pass through vcf without intersect"}

outputs:
  intersected_vcf:
    type: 'File?'
    outputBinding:
      glob: '*.bed_intersect.vcf.gz'
    secondaryFiles: ['.tbi']
