cwlVersion: v1.0
class: CommandLineTool
id: gatk4_getpileupsummary
label: GATK Pileup
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: ${ return inputs.max_memory * 1000 }
    coresMin: 2
baseCommand: [/gatk, GetPileupSummaries]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m"
      -I $(inputs.aligned_reads.path)
      -V $(inputs.exac_common_vcf.path)
      -L $(inputs.interval_list.path)
      -R $(inputs.reference.path)
      -O $(inputs.aligned_reads.nameroot).pileupsummary.table

inputs:
  aligned_reads: {type: 'File', secondaryFiles: ['.crai']}
  reference: File
  interval_list: File
  exac_common_vcf: {type: 'File', secondaryFiles: [.tbi]}
  max_memory: {type: 'int?', default: 2, doc: "Maximum memory in GB for GATK GetPileupSummaries to use"}

outputs:
  pileup_table:
    type: File
    outputBinding:
      glob: '*.pileupsummary.table'

