cwlVersion: v1.0
class: CommandLineTool
id: samtools_calmd
doc: "Recalculate MD tags from input"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 12000
    coresMin: $(inputs.threads)

baseCommand: ["/bin/bash -c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      samtools calmd -@ $(inputs.threads) -b --reference $(inputs.reference.path) $(inputs.input_reads.path) > $(inputs.input_reads.nameroot).bam;
      samtools index -@ $(inputs.threads) $(inputs.input_reads.nameroot).bam $(inputs.input_reads.nameroot).bai

inputs:
  input_reads: File
  threads:
    type: ['null', int]
    default: 16
  reference: {type: File, secondaryFiles: [.fai]}
outputs:
  bam_file:
    type: File
    outputBinding:
      glob: '*.bam'
    secondaryFiles: [^.bai]
