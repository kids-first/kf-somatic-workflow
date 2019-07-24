cwlVersion: v1.0
class: CommandLineTool
id: cnvkit-batch
requirements: 
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement 
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/milos_nikolic/cnvkit:0.9.3'
  - class: ResourceRequirement
    ramMin: 32000
    coresMin: 16
arguments: 
  - position: 1
    shellQuote: false
    valueFrom: >-
      ln -s $(inputs.normal_bam.path) .; ln -s $(inputs.normal_bam.secondaryFiles[0].path) ./$(inputs.normal_bam.basename).bai

      ln -s $(inputs.tumor_bam.path) .; ln -s $(inputs.tumor_bam.secondaryFiles[0].path) ./$(inputs.tumor_bam.basename).bai

      cnvkit.py batch $(inputs.tumor_bam.basename)
      --normal $(inputs.normal_bam.basename)
      --fasta $(inputs.reference_fasta.path)
      --targets $(inputs.target_regions.path)
      --output-reference $(inputs.output_reference_name) 
      --diagram --scatter

      cnvkit.py call $(inputs.tumor_bam.nameroot).cns -o $(inputs.tumor_bam.nameroot).call.cns

      cnvkit.py export vcf $(inputs.tumor_bam.nameroot).call.cns -o $(inputs.tumor_bam.nameroot).vcf

      mv $(inputs.tumor_bam.basename).call.cns $(inputs.output_basename).calls.cns

      mv $(inputs.tumor_bam.basename).diagram.pdf $(inputs.output_basename).diagram.pdf

      mv $(inputs.tumor_bam.basename).scatter.pdf $(inputs.output_basename).scatter.pdf

      mv $(inputs.tumor_bam.basename).vcf $(inputs.output_basename).cnvkit.vcf
      
      gzip $(inputs.output_basename).cnvkit.vcf

inputs:
  tumor_bam: {type: File, doc: "tumor bam file", secondaryFiles: [.bai]}
  normal_bam: {type: File, doc: "matched normal bam file", secondaryFiles: [.bai]}
  reference_fasta: {type: File, doc: "reference fasta file", secondaryFiles: [.fai]}
  target_regions: {type: File, doc: "target intervals"}
  output_reference_name: string
  output_basename: string

outputs:
  output_vcf:
    type: File
    outputBinding:
      glob: '*.vcf'
    secondaryFiles: [.tbi]
  output_calls:
    type: File
    outputBinding:
      glob: '*.call.cns'
  output_scatter:
    type: File
    outputBinding:
      glob: '*.scatter.pdf'
  output_diagram:
    type: File
    outputBinding:
      glob: '*.diagram.pdf'
