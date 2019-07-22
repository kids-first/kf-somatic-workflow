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
      cnvkit.py batch $(inputs.tumor_bam.path)
      --normal $(inputs.normal_bam.path)
      --fasta $(inputs.reference_fasta.path)
      --targets $(inputs.target_regions.path)
      --output-reference $(inputs.output_reference_name)
      --diagram
      --scatter

      cnvkit.py call $(inputs.tumor_bam.basename) -o $(inputs.tumor_bam.basename).call.cns

      cnvkit.py export vcf $(inputs.tumor_bam.basename).call.cns -o $(tumor_bam.basename).vcf

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
      glob: '*.CNV.vcf.gz'
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
      glob: '*diagram.pdf'
