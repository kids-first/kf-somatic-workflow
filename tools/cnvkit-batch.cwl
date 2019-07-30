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
      ln -s $(inputs.tumor_bam.path) .; ln -s $(inputs.tumor_bam.secondaryFiles[0].path) ./$(inputs.tumor_bam.basename).bai
      ln -s $(inputs.normal_bam.path) .; ln -s $(inputs.normal_bam.secondaryFiles[0].path) ./$(inputs.normal_bam.basename).bai
      
      cnvkit.py batch $(inputs.tumor_bam.path) 
      --normal $(inputs.normal_bam.path)
      --fasta $(inputs.reference.path) 
      --targets $(inputs.capture_regions.path) 
      --annotate $(inputs.annotation_file.path)
      --output-reference $(inputs.output_basename)_cnvkit_reference.cnn 
      --diagram 
      --scatter

      mv $(inputs.tumor_bam.nameroot)-diagram.pdf $(inputs.output_basename).diagram.pdf
      mv $(inputs.tumor_bam.nameroot)-scatter.pdf $(inputs.output_basename).scatter.pdf

      cnvkit.py call $(inputs.tumor_bam.nameroot).cns 
      -o $(inputs.output_basename).call.cns
      
      cnvkit.py export vcf $(inputs.output_basename).call.cns -o $(inputs.output_basename).vcf

      cnvkit.py metrics $(inputs.output_basename).cnr -s $(inputs.tumor_bam.nameroot).cns
      -o $(inputs.output_basename).metrics.txt

      cnvkit.py gainloss $(inputs.output_basename).cnr -o $(inputs.output_basename).gainloss.txt


inputs:
  tumor_bam: {type: File, doc: "tumor bam file", secondaryFiles: [.bai]}
  normal_bam: {type: File, doc: "normal bam file", secondaryFiles: [.bai]}
  reference: {type: File, doc: "fasta file", secondaryFiles: [.fai]}
  capture_regions: {type: File, doc: "target regions"}
  annotation_file: {type: File, doc: "refFlat.txt file"}
  output_basename: string

outputs:
  output_cnr: 
    type: File
    outputBinding:
      glob: '*.cnr'
  output_vcf:
    type: File
    outputBinding:
      glob: '*.vcf'
  output_calls:
    type: File
    outputBinding:
      glob: '*.call.cns'
  output_scatter:
    type: File
    outputBinding:
      glob: '*scatter.pdf'
  output_diagram:
    type: File
    outputBinding:
      glob: '*diagram.pdf'
  output_metrics: 
    type: File
    outputBinding:
      glob: '*.metrics.txt'
  output_gainloss: 
    type: File
    outputBinding:
      glob: '*.gainloss.txt'
