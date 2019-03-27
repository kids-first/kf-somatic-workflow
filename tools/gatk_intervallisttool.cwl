cwlVersion: v1.0
class: CommandLineTool
id: gatk4_Mutect2
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/gatk:4.0.12.0'
  - class: ResourceRequirement
    ramMin: 2000
baseCommand: [/gatk, IntervalListTools]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      --java-options "-Xmx2000m"
      --SCATTER_COUNT=50
      --SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW
      --UNIQUE=true --SORT=true --BREAK_BANDS_AT_MULTIPLES_OF=1000000
      --INPUT=$(inputs.interval_list.path) --OUTPUT=.
inputs:
  interval_list: File
outputs:
  output:
    type: 'File[]'
    outputBinding:
      glob: 'temp*/*.interval_list'
