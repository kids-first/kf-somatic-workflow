cwlVersion: v1.0
class: CommandLineTool
id: kf-vep-somatic-annotate
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 24000
  - class: DockerRequirement
    dockerPull: 'migbro/vep:r93'
baseCommand: [tar, -xzfv ]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      perl /ensembl-vep/vep
      --cache --dir_cache $PWD
      --cache_version 93
      --vcf
      --symbol
      --canonical
      --variant_class
      --offline
      --hgvs
      --hgvsg
      --fork 7
      --sift b
      -i $(inputs.input_vcf.path)
      -o $(inputs.output_basename).vep.vcf
      --fasta Homo_sapiens_assembly38.fasta

inputs:
  reference: { type: File,  secondaryFiles: [.fai] }
  input_vcf:
    type: File
    secondaryFiles: [tbi]
  output_basename: string
  cache: {type: File }

outputs:
  output:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: [.tbi]