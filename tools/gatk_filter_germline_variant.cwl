cwlVersion: v1.0
class: CommandLineTool
id: kfdrc-gatk_variantfiltration
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 2
    coresMax: 4
  - class: DockerRequirement
    dockerPull: 'kfdrc/gatk:4.1.1.0'
baseCommand: [/gatk]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      --java-options "-Xmx7000m"
      SelectVariants
      --exclude-filtered TRUE
      -select-type SNP
      -V $(inputs.input_vcf.path)
      -O snp_pass.vcf.gz

      /gatk VariantFiltration
      --java-options "-Xmx7000m"
      --filter-name GATK_QD
      --filter-expression "QD < 2.0"
      --filter-name GATK_FS
      --filter-expression "FS > 60.0"
      --filter-name GATK_MQ
      --filter-expression "MQ < 40.0"
      --filter-name GATK_MQRankSum
      --filter-expression "MQRankSum < -12.5"
      --filter-name GATK_ReadPosRankSum
      --filter-expression "ReadPosRankSum < -8.0"
      --filter-name KFDRC_DP10
      --filter-expression "DP < 10"
      -V snp_pass.vcf.gz -O $(inputs.output_basename).gatk.hardfiltered.vcf.gz

      /gatk
      SelectVariants
      --exclude-filtered TRUE
      -V $(inputs.output_basename).gatk.hardfiltered.vcf.gz
      -O $(inputs.output_basename).gatk.hardfiltered.PASS.vcf.gz

inputs:
  reference_fasta: File
  input_vcf: {type: File, secondaryFiles: [.tbi]}
  output_basename: string
outputs:
  filtered_vcf:
    type: File
    outputBinding:
      glob: '*.gatk.hardfiltered.vcf.gz'
    secondaryFiles: [.tbi]
  filtered_pass_vcf:
    type: File
    outputBinding:
      glob: '*.gatk.hardfiltered.PASS.vcf.gz'
    secondaryFiles: [.tbi]

