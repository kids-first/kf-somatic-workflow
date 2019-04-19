cwlVersion: v1.0
class: CommandLineTool
id: gatk4_intervallist2bed
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: 2000

baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail

      /gatk IntervalListTools
      --java-options "-Xmx2000m"
      --SCATTER_COUNT=200
      --SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW
      --UNIQUE=true --SORT=true --BREAK_BANDS_AT_MULTIPLES_OF=1000000
      --INPUT=$(inputs.interval_list.path) --OUTPUT=.
      && seq -w 200 | xargs -I N -P 4
      /gatk IntervalListToBed --java-options -Xmx100m -I temp_0N_of_200/scattered.interval_list -O temp_0N_of_200/scattered.interval_list.N.bed
inputs:
  interval_list: File
outputs:
  output:
    type: 'File[]'
    outputBinding:
      glob: 'temp*/*.bed'
