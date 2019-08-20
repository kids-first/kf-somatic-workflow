## this tool takes in panel of normal. must run the cnvkit reference tool on cavatica first
cwlVersion: v1.0
class: CommandLineTool
id: cnvkit-calling-pipeline
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
      
      cnvkit.py coverage $(inputs.tumor_bam.path) 
      $(inputs.target_coverage.path) 
      -o $(inputs.output_basename).targetcoverage.cnn
      
      cnvkit.py coverage $(inputs.tumor_bam.path) 
      $(inputs.antitarget_coverage.path) 
      -o $(inputs.output_basename).antitargetcoverage.cnn
     
      cnvkit.py fix $(inputs.output_basename).targetcoverage.cnn 
      $(inputs.output_basename).antitargetcoverage.cnn 
      $(inputs.reference_coverage_file.path) 
      -o $(inputs.output_basename).cnv.cnr

      cnvkit.py segment $(inputs.output_basename).cnv.cnr 
      --drop-low-coverage 
      -o $(inputs.output_basename).cns

      cnvkit.py export vcf $(inputs.output_basename).cns -o $(inputs.output_basename).vcf

      cnvkit.py call $(inputs.output_basename).cns 
      -v $(inputs.output_basename).vcf
      -o $(inputs.output_basename).call.cns

      cnvkit.py scatter $(inputs.output_basename).cnv.cnr
      --segment $(inputs.output_basename).cns
      --y-max $(inputs.y_max)
      --y-min $(inputs.y_min)
      -o $(inputs.output_basename).scatter.pdf

      cnvkit.py diagram $(inputs.output_basename).cnv.cnr
      --segment $(inputs.output_basename).cns
      -o $(inputs.output_basename).diagram.pdf


inputs:
  tumor_bam: {type: File, doc: "tumor bam file", secondaryFiles: [.bai]}
  target_coverage: {type: File, doc: "targetcoverage.bed from CNVkit reference wf"}
  antitarget_coverage: {type: File, doc: "antitargetcoverage.bed from CNVkit reference wf"}
  reference_coverage_file: {type: File, doc: "reference file from CNVkit reference wf"}
  output_basename: string
  y_min: {type: ['null', string], doc: "y min for scatter plot"}
  y_max: {type: ['null', string], doc: "y max for scatter plot"}

outputs:
  output_cnr: 
    type: File
    outputBinding:
    glob: '*.cnv.cnr'
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
      glob: '*.scatter.pdf'
  output_diagram:
    type: File
    outputBinding:
      glob: '*.diagram.pdf'
