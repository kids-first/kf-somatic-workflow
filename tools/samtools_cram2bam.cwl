cwlVersion: v1.0
class: CommandLineTool
id: samtools_cram2bam
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 12000
    coresMin: $(inputs.threads)
  
baseCommand: [samtools, view]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -@ $(inputs.threads)
      -bh $(inputs.input_reads.path)
      -T $(inputs.reference.path)
      > $(inputs.input_reads.nameroot).bam
      && samtools index $(inputs.input_reads.nameroot).bam $(inputs.input_reads.nameroot).bai
inputs:
  input_reads: File
  threads:
    type: ['null', int]
    default: 16
  reference: File
outputs:
  bam_file:
    type: File
    outputBinding:
      glob: '*.bam'
    secondaryFiles: [^.bai]
