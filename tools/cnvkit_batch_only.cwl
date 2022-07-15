cwlVersion: v1.2
class: CommandLineTool
id: cnvkit-batch-only
doc: "This tools is a bit more restrictive than the cnvkit_batch.cwl, meant to be run to get needed .cns file when needed"
requirements: 
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement 
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/milos_nikolic/cnvkit:0.9.3'
  - class: ResourceRequirement
    ramMin: 32000
    coresMin: $(inputs.threads)
  - class: InitialWorkDirRequirement
    listing: [$(inputs.input_sample),$(inputs.input_control)]
arguments: 
  - position: 1
    shellQuote: false
    valueFrom: >-
      cnvkit.py batch
      ${
        if (inputs.cnvkit_cnn == null){
          return "--output-reference " + inputs.output_basename + "_cnvkit_reference.cnn";
        }
        else{
            return "";
        }
      }
inputs:
  input_sample: { type: File, doc: "tumor bam file", secondaryFiles: [^.bai], inputBinding: { position: 1}}
  input_control: { type: 'File?', doc: "normal bam file - can skip in .cnn file supplied", secondaryFiles: ['^.bai?', '.crai?'], inputBinding: { position: 2, prefix: "--normal"} }
  reference: { type: 'File?', doc: "fasta file, needed if cnv kit cnn not already built", secondaryFiles: [.fai], inputBinding: { position: 2, prefix: "--fasta"} }
  cnvkit_cnn: { type: 'File?', doc: "If running using an existing .cnn, supply here", inputBinding: { position: 2, prefix: "--reference" } }
  capture_regions: { type: 'File?', doc: "target regions for WES", inputBinding: { prefix: '--targets', position: 2 } }
  annotation_file: { type: 'File?', doc: "refFlat.txt file, needed if cnv kit cnn not already built", inputBinding: { position: 2, prefix: "--annotate"} }
  output_basename: string
  run_mode: { type: ['null', {type: enum, name: wgs_mode, symbols: ["wgs", "hybrid", "amplicon"]}], doc: "Choose rum method", default: "wgs", inputBinding: { position: 2, prefix: "-m"} }
  threads: { type: 'int?', doc: 'Num threads to use', default: 16, inputBinding: { position: 2, prefix: "-p"} }
  scatter_plot: { type: 'boolean?', inputBinding: { prefix: "--scatter", position: 2 }, doc: "Create a whole-genome copy ratio profile as a PDF scatter plot." }
  diagram_plot: { type: 'boolean?', inputBinding: { prefix: "--diagram", position: 2 }, doc: "Create an ideogram of copy ratios on chromosomes as a PDF." }
  male_input_flag: { type: 'boolean?', doc: "Is the input male?", default: false, inputBinding: { position: 2,  prefix: "--male-reference"} }
outputs:
  output_cnr: 
    type: File
    outputBinding:
      glob: '*.cnr'
  output_cns:
    type: 'File'
    outputBinding:
      glob: '*.cns'
  output_scatter:
    type: 'File?'
    outputBinding:
      glob: '*-scatter.pdf'
  output_diagram:
    type: 'File?'
    outputBinding:
      glob: '*-diagram.pdf'
  output_cnn:
    type: 'File?'
    outputBinding:
      glob: '*_reference.cnn'
    doc: "Output if starting from cnn scratch. Should not appear if an existing .cnn was given as input."

