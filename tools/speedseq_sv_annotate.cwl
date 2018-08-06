cwlVersion: v1.0
class: CommandLineTool
id: speedseq_sv_annotate
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 8000
  - class: DockerRequirement
    dockerPull: 'speedseq:latest'
baseCommand: [/speedseq/bin/speedseq, sv]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -dB $(inputs.input_align.path)
      -R $(inputs.reference.path)
      -t 8
      -o $(inputs.output_basename)
inputs:
  reference: { type: File,  secondaryFiles: [.fai] }
  input_align: { type: File,  secondaryFiles: [.crai|.cram] }
  output_basename: string
outputs:
  output:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: [.tbi]