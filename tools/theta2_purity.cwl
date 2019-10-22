cwlVersion: v1.0
class: CommandLineTool
id: theta2
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/theta2:0.7'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 32000
    coresMin: 8

baseCommand: [/THetA/bin/RunTHetA]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      --TUMOR_FILE $(inputs.tumor_snp.path)
      --NORMAL_FILE $(inputs.normal_snp.path)
      --OUTPUT_PREFIX $(inputs.output_basename)
      --NUM_PROCESSES 8
      --MIN_FRAC $(inputs.min_frac)
      $(inputs.interval_count.path)
inputs:
  tumor_snp: File
  normal_snp: File
  interval_count: File
  output_basename: string
  min_frac: {type: ['null', float], doc: "Minimum fraction of genome with copy umber alterations.  Default is 0.05", default: 0.05}
outputs:
  n2_graph:
    type: File
    outputBinding:
      glob: '*.n2.graph.pdf'
  n2_results:
    type: File
    outputBinding:
      glob: '*.n2.results'
  n2_withBounds:
    type: File
    outputBinding:
      glob: '*.n2.withBounds'
  n3_graph:
    type: File
    outputBinding:
      glob: '*.n3.graph.pdf'
  n3_results:
    type: File
    outputBinding:
      glob: '*.n3.results'
  n3_withBounds:
    type: File
    outputBinding:
      glob: '*.n3.withBounds'
  best_results:
    type: File
    outputBinding:
      glob: '*.BEST.results'
