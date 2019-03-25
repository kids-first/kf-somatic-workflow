cwlVersion: v1.0
class: CommandLineTool
id: seurat
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'migbro/seurat:2.6'
  - class: ResourceRequirement
    ramMin: 3000
baseCommand: [java]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      "-Xms3000m -Xmx3000m"
      -jar
      -R $(inputs.reference.path)
      -I:dna_tumor $(inputs.input_tumor_bam.path)
      -I:dna_normal $(inputs.input_normal_bam.path)
      --indels
      -go $(inputs.input_tumor_bam.nameroot).$(inputs.interval_list.nameroot).large_events.txt
      -L $(inputs.interval_list.path)
      -refseq $(inputs.refseq_rod.path)
      -o $(inputs.input_tumor_bam.nameroot).$(inputs.interval_list.nameroot).somatic.seurat.vcf

inputs:
  reference: {type: File, secondaryFiles: [^.dict, .fai]}
  refseq_rod: type: File
  input_tumor_bam: {type: File, secondaryFiles: [^.bai]}
  input_tumor_name: string
  input_normal_bam: {type: File, secondaryFiles: [^.bai]}
  input_normal_name: string
  interval_list: File

outputs:
  seurat_vcf:
    type: File
    outputBinding:
      glob: '*.vcf'
  seurat_events:
    type: File
    outputBinding:
      glob: '*.txt'