cwlVersion: v1.0
class: CommandLineTool
id: strelka2
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 10000
    coresMin: 36
  - class: DockerRequirement
    dockerPull: 'obenauflab/strelka'

baseCommand: [/strelka-2.9.3.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      --normalBam $(inputs.input_normal_cram.path)
      --tumorBam $(inputs.input_tumor_cram.path)
      --ref $(inputs.reference.path)
      --callRegions $(inputs.hg38_strelka_bed.path)
      --runDir=./ && ./runWorkflow.py
      -m local
      -j 36

inputs:
    reference: { type: File, secondaryFiles: [^.dict, .fai] }
    hg38_strelka_bed: { type: File, secondaryFiles: [.tbi], label: bed file of hg38 chrs without contigs }
    input_tumor_cram: { type: File, secondaryFiles: [.crai] }
    input_normal_cram: { type: File, secondaryFiles: [.crai] }
outputs:
  - id: output_snv
    type: File
    outputBinding:
      glob: 'results/variants/*.snvs.vcf.gz'
    secondaryFiles: [.tbi]
  - id: output_indel
    type: File
    outputBinding:
      glob: 'results/variants/*.indels.vcf.gz'
    secondaryFiles: [.tbi]

