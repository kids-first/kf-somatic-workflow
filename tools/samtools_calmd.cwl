cwlVersion: v1.0
class: CommandLineTool
id: samtools_calmd
doc: "Recalculate MD tags from input"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/brownm28/samtools:1.20-multi-arch'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 12000
    coresMin: $(inputs.threads)

baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      . /link_samtools.sh;
      samtools calmd -@ $(inputs.threads) -b --reference $(inputs.reference.path) $(inputs.input_reads.path) > $(inputs.input_reads.nameroot).bam;
      samtools index -@ $(inputs.threads) $(inputs.input_reads.nameroot).bam $(inputs.input_reads.nameroot).bai

inputs:
  input_reads: File
  threads:
    type: ['null', int]
    default: 8
  reference: {type: File, secondaryFiles: [.fai]}
outputs:
  bam_file:
    type: File
    outputBinding:
      glob: '*.bam'
    secondaryFiles: [^.bai]
