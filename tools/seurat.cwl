cwlVersion: v1.0
class: CommandLineTool
id: seurat
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/seurat:latest'
  - class: ResourceRequirement
    ramMin: 3000
    coresMin: 4
baseCommand: [java]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -Xmx3000m
      -jar
      /Seurat-2.6.jar
      -T Seurat
      --filter_mismatching_base_and_quals
      -R $(inputs.reference.path)
      -I:dna_tumor $(inputs.input_tumor_bam.path)
      -I:dna_normal $(inputs.input_normal_bam.path)
      --indels
      -L $(inputs.interval_list.path)
      -o $(inputs.input_tumor_bam.nameroot).$(inputs.interval_list.nameroot).somatic.seurat.vcf
      && bgzip $(inputs.input_tumor_bam.nameroot).$(inputs.interval_list.nameroot).somatic.seurat.vcf
      && tabix $(inputs.input_tumor_bam.nameroot).$(inputs.interval_list.nameroot).somatic.seurat.vcf.gz

inputs:
  reference: {type: File, secondaryFiles: [^.dict, .fai]}
  input_tumor_bam: {type: File, secondaryFiles: [^.bai]}
  input_normal_bam: {type: File, secondaryFiles: [^.bai]}
  interval_list: File

outputs:
  seurat_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: [.tbi]
