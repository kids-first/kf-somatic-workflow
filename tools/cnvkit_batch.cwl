## this tool is made to run wgs or wes, germline/somatic or paired, using the batch method (no panel of normals needed)
cwlVersion: v1.2
class: CommandLineTool
id: cnvkit-batch
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
      touch $(inputs.input_sample.secondaryFiles[0].path)

      $(inputs.input_control != null ? "touch " + inputs.input_control.secondaryFiles[0].path : "")

      cnvkit.py batch
      --diagram 
      --scatter
      ${
        if (inputs.cnvkit_cnn == null){
          var arg = "--output-reference " + inputs.output_basename + "_cnvkit_reference.cnn --fasta " + inputs.reference.path + " --annotate " + inputs.annotation_file.path;
          if (inputs.input_control != null) {
              arg += " --normal " + inputs.input_control.path;
          }
        }
        else{
          var arg = "--reference " + inputs.cnvkit_cnn.path;
          var msex = ['m','y','male','Male']
          if (msex.indexOf(inputs.sex) >= 0){
            arg += " --male-reference";
          }
        }
        return arg;
      }
  - position: 10
    shellQuote: false
    valueFrom: >-
      && cnvkit.py call
      $(inputs.input_sample.nameroot).cns
      ${
        var arg = "--sample-sex " + inputs.sex;
        var msex = ['m','y','male','Male']
        if (msex.indexOf(inputs.sex) >= 0){
          arg += " --male-reference";
        }
        return arg;
      }
      -o $(inputs.output_basename).call.cns
  - position: 20
    shellQuote: false
    valueFrom: >-
      && cnvkit.py export seg $(inputs.output_basename).call.cns | sed -e 's/$(inputs.output_basename).call/$(inputs.tumor_sample_name)/' > $(inputs.output_basename).call.seg
      && cnvkit.py metrics $(inputs.input_sample.nameroot).cnr -s $(inputs.input_sample.nameroot).cns
      -o $(inputs.output_basename).metrics.txt
      && cnvkit.py gainloss $(inputs.input_sample.nameroot).cnr -o $(inputs.output_basename).gainloss.txt
      && mv $(inputs.input_sample.nameroot).cnr $(inputs.output_basename).cnr
      && mv $(inputs.input_sample.nameroot)-diagram.pdf $(inputs.output_basename).diagram.pdf
      && mv $(inputs.input_sample.nameroot)-scatter.pdf $(inputs.output_basename).scatter.pdf

inputs:
  # batch params
  threads: {type: 'int?', default: 16, inputBinding: { position: 2, prefix: "-p" } }
  wgs_mode: { type: ['null', {type: enum, name: wgs_mode, symbols: ["hybrid", "wgs", "amplicon"]}],
    default: "hybrid",
    doc: "for WGS mode, input wgs. leave blank for hybrid mode",
    inputBinding: { position: 2, prefix: "-m"} }
  access: { type: 'File?', doc: "Regions of accessible sequence on chromosomes (.bed), as output by the 'access' command.",
    inputBinding: { position: 2, prefix: "--access"} }
  input_sample: {type: File, doc: "tumor bam file", secondaryFiles: [^.bai],
    inputBinding: { position: 3 } }
  capture_regions: {type: ['null', File], doc: "target regions for WES",
    inputBinding: {position: 2, prefix: "--targets"} }
  # call
  b_allele_vcf: {type: ['null', File], doc: "b allele germline vcf, if available",
    inputBinding: { position: 11, prefix: "--vcf"} }
  input_control: {type: ['null', File], doc: "normal bam file", secondaryFiles: [^.bai]}
  # custom
  reference: {type: ['null', File], doc: "fasta file, needed if cnv kit cnn not already built", secondaryFiles: [.fai]}
  cnvkit_cnn: {type: ['null', File], doc: "If running using an existing .cnn, supply here"}
  annotation_file: {type: ['null', File], doc: "refFlat.txt file, needed if cnv kit cnn not already built"}
  output_basename: string
  tumor_sample_name: string
  sex: {type: ['null', {type: enum, name: sex, symbols: ["x", "y"]}], doc: "Sex, for simplicity x for female y for male", default: "x"}
outputs:
  output_cnr: 
    type: File
    outputBinding:
      glob: '*.cnr'
  output_cns:
    type: File
    outputBinding:
      glob: '$(inputs.input_sample.nameroot).cns'
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
