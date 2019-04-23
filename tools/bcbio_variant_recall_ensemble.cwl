cwlVersion: v1.0
class: CommandLineTool
id: bcbio_vr
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 6000
    coresMin: 8
  - class: DockerRequirement
    dockerPull: 'kfdrc/bcbio_vr:0.2.4'

baseCommand: [/bcbio-variation-recall]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      ensemble --cores 8 --numpass $(inputs.min_overlap)
      --names $(inputs.tool_name_csv) $(inputs.output_basename).caller_consensus.vcf.gz
      $(inputs.reference.path)
inputs:
    reference: File
    tool_name_csv: string
    output_basename: string
    min_overlap: {type: int, doc: "Min number of callers to declare consensus.  Default is 2"}
    input_vcfs:
      type:
        type: array
        items: File
      secondaryFiles: [.tbi]
      inputBinding:
        position: 1

outputs:
  consensus_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: ['.tbi']
