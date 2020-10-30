cwlVersion: v1.0
class: CommandLineTool
id: samtools_reheader_bam
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9'
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4
  - class: InlineJavascriptRequirement
baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -euxo pipefail

      samtools view -H $(inputs.input_reads.path) | sed -E "s/SM:\S+/SM:$(inputs.sample_name)/" | samtools reheader - $(inputs.input_reads.path) > $(inputs.input_reads.basename)

      samtools index -@ 4 $(inputs.input_reads.basename) $(inputs.input_reads.nameroot).bai

inputs:
  input_reads: {type: File, secondaryFiles: [^.bai]}
  sample_name: string
outputs:
  bam_file:
    type: File
    outputBinding:
      glob: '*.bam'
    secondaryFiles: [^.bai]
