cwlVersion: v1.2
class: CommandLineTool
id: samtools_calmd
doc: "Recalculate MD tags from input"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'staphb/samtools:1.21'
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.cpu)

baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      samtools calmd -@ $(inputs.cpu) -b --reference $(inputs.reference.path) $(inputs.input_reads.path) > $(inputs.input_reads.nameroot).bam;
      samtools index -@ $(inputs.cpu) $(inputs.input_reads.nameroot).bam $(inputs.input_reads.nameroot).bai

inputs:
  input_reads: { type: 'File', secondaryFiles: [{pattern: '^.bai', required: false}, {pattern: '.bai', required: false}, {pattern: '.crai', required: false}, {pattern: '^.crai', required: false}]}
  reference: {type: 'File', secondaryFiles: [{pattern: '.fai', required: true}]}
  cpu: { type: 'int?', default: 8, doc: "Number of CPUs for this task." }
  ram: { type: 'int?', default: 16, doc: "GB of RAM for this task." }
outputs:
  bam_file:
    type: File
    outputBinding:
      glob: '*.bam'
    secondaryFiles: [{pattern: '^.bai', required: true}]
