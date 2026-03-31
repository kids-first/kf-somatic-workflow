cwlVersion: v1.2
class: CommandLineTool
id: bcftools-gen-b-allele
doc: "Applies GATK germline hard filters rules for b allele output"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 16000
    coresMin: $(inputs.threads)
  - class: DockerRequirement
    dockerPull: 'staphb/bcftools:1.20'

baseCommand: [bcftools, view]
arguments:
  - position: 2
    shellQuote: false
    valueFrom: >-
      -v snps
      -f "PASS,."
      -e "INFO/QD < 2.0 ||
      INFO/FS > 60.0 ||
      INFO/MQ < 40.0 ||
      INFO/MQRankSum < -12.5 ||
      INFO/ReadPosRankSum < -8.0 ||
      INFO/DP < 10"
      -O z
      -o $(inputs.output_basename).$(inputs.tool_name).hardfiltered.PASS.vcf.gz
      --write-index=tbi

inputs:
    input_vcf: {type: 'File', secondaryFiles: ['.tbi'],
      inputBinding: { position: 3} }
    threads: {type: 'int?', doc: "Number of compression/decompression threads",
      default: 8, inputBinding: {position: 1, prefix: "--threads"} }
    output_basename: { type: 'string?', default: "b_allele_cnv" }
    tool_name: { type: 'string?', default: "gatk" }

outputs:
  bcftools_b_allele_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: ['.tbi']
