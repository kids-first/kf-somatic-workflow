cwlVersion: v1.0
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
          var arg = "--output-reference " + inputs.output_basename + "_cnvkit_reference.cnn --fasta " + inputs.reference.path + " --annotate " + inputs.annotation_file.path;
          if (inputs.input_control != null) {
              arg += " --normal " + inputs.input_control.path;
          }
        }
        else{
          var arg = "--reference " + inputs.cnvkit_cnn.path;

        }
        return arg;
      }
      --diagram
      --scatter

inputs:
  input_sample: {type: File, doc: "tumor bam file", secondaryFiles: [^.bai], inputBinding: { position 1}}
  input_control: {type: ['null', File], doc: "normal bam file", secondaryFiles: [^.bai]}
  reference: {type: ['null', File], doc: "fasta file, needed if cnv kit cnn not already built", secondaryFiles: [.fai]}
  cnvkit_cnn: {type: ['null', File], doc: "If running using an existing .cnn, supply here"}
  capture_regions: {type: ['null', File], doc: "target regions for WES", inputBinding: { prefix: '--targets', position: 2 }}
  annotation_file: {type: ['null', File], doc: "refFlat.txt file, needed if cnv kit cnn not already built"}
  output_basename: string
  wgs_mode: { type: 'boolean?', doc: "for WGS mode, set to true.", default: true, inputBinding: { position: 2, prefix: '-m wgs', shellQuote: false } }
  threads: { type: 'int?', doc: 'Num threads to use', default: 16, inputBinding: { position: 2, prefix: "-p"} }
  sex: {type: 'boolean?', doc: "Is the input male?", default: false, inputBinding: { position: 2,  prefix: "--male-reference"}}
outputs:
  output_cnr: 
    type: File
    outputBinding:
      glob: '*.cnr'
  output_cns:
    type: 'File[]'
    outputBinding:
      glob: '*.cns'
  output_scatter:
    type: File
    outputBinding:
      glob: '*.scatter.pdf'
  output_diagram:
    type: File
    outputBinding:
      glob: '*.diagram.pdf'
  output_metrics: 
    type: File
    outputBinding:
      glob: '*.metrics.txt'
  output_gainloss: 
    type: File
    outputBinding:
      glob: '*.gainloss.txt'
  output_cnn:
    type: ['null', File]
    outputBinding:
      glob: '*_reference.cnn'
    doc: "Output if starting from cnn scratch.  Should not appear if an existing .cnn was given as input."
  output_seg:
    type: File
    outputBinding:
      glob: '*.seg'
