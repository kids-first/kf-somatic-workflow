cwlVersion: v1.0
class: CommandLineTool
id: kfdrc-vep-somatic-annotate
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 24000
    coresMin: 16
  - class: DockerRequirement
    dockerPull: 'kfdrc/vep:r93'
baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail

      tar -xzf $(inputs.cache.path)
      && perl /ensembl-vep/vep
      --cache --dir_cache $PWD
      --cache_version 93
      --vcf
      --symbol
      --canonical
      --variant_class
      --offline
      --hgvs
      --hgvsg
      --fork 15
      --sift b
      --vcf_info_field ANN
      -i $(inputs.input_vcf.path)
      -o STDOUT
      --stats_file $(inputs.output_basename)_stats.$(inputs.tool_name).html
      --warning_file $(inputs.output_basename)_warnings.$(inputs.tool_name).txt
      --fasta $(inputs.reference.path) |
      /ensembl-vep/htslib/bgzip -c > $(inputs.output_basename).$(inputs.tool_name).vep.vcf.gz
      && /ensembl-vep/htslib/tabix $(inputs.output_basename).$(inputs.tool_name).vep.vcf.gz

inputs:
  reference: { type: File,  secondaryFiles: [.fai], label: Fasta genome assembly with index }
  input_vcf:
    type: File
    secondaryFiles: [.tbi]
  output_basename: string
  tool_name: string
  cache: { type: File, label: tar gzipped cache from ensembl/local converted cache }

outputs:
  output_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
  output_tbi:
    type: File
    outputBinding:
      glob: '*.vcf.gz.tbi'

  output_html:
    type: File
    outputBinding:
      glob: '*.html'
  warn_txt:
    type: ["null", File]
    outputBinding:
      glob: '*.txt'
