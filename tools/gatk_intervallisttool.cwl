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
      --SCATTER_COUNT=$(inputs.scatter_ct)
      --SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW
      --UNIQUE=true --SORT=true --BREAK_BANDS_AT_MULTIPLES_OF=$(inputs.bands)
      --INPUT=$(inputs.interval_list.path) --OUTPUT=.
      && CT=`find . -name 'temp_0*' | wc -l`
      && seq -f "%04g" $CT | xargs -I N -P 4
      /gatk IntervalListToBed --java-options -Xmx100m -I temp_N_of_$CT/scattered.interval_list -O temp_N_of_$CT/scattered.interval_list.N.bed
inputs:
  interval_list: File
  bands: int
  scatter_ct: int
outputs:
  output:
    type: 'File[]'
    outputBinding:
      glob: 'temp*/*.bed'
