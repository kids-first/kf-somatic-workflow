## this tool is made to run wgs or wes, germline/somatic or paired, using the batch method (no panel of normals needed)
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
    coresMin: $(inputs.threads)
arguments: 
  - position: 1
    shellQuote: false
    valueFrom: >-
      ln -s $(inputs.input_sample.path) .; ln -s $(inputs.input_sample.secondaryFiles[0].path) ./$(inputs.input_sample.basename).bai
      
      ${ 
          var cmd = "";
          if (inputs.input_control != null) {
              cmd = "ln -s " + inputs.input_control.path + " .; ln -s " + inputs.input_control.secondaryFiles[0].path + " ./" + inputs.input_control.basename + ".bai"
          }
          return cmd;
      }

      cnvkit.py batch
      -p $(inputs.threads)
      ${
          var cmd = "";
          if (inputs.wgs_mode == 'Y') {
              cmd = " -m wgs ";
          }
          return cmd;
      }
      $(inputs.input_sample.path) 
      ${
          var cmd = "--normal ";
          if (inputs.input_control != null) {
              cmd += inputs.input_control.path + " ";
          }
          return cmd;
      }
      --fasta $(inputs.reference.path) 
      ${
          var cmd = "";
          if (inputs.capture_regions != null) {
              cmd = "--targets " + inputs.capture_regions.path;
          }
          return cmd;
      }
      --drop-low-coverage
      --annotate $(inputs.annotation_file.path)
      ${
        var arg = "--output-reference " + inputs.output_basename + "_cnvkit_reference.cnn";
        if(inputs.cnv_kit_cnn != null){
          arg = "--reference " + inputs.cnv_kit_cnn.path;
        }
        return arg;
      }
      --diagram 
      --scatter

      cnvkit.py call $(inputs.input_sample.nameroot).cns
      ${
        var arg = "";
        if (inputs.b_allele_vcf != null){
          arg = "--vcf " + inputs.b_allele_vcf.path;
        }
        return arg;
      }
      -o $(inputs.output_basename).call.cns
      
      cnvkit.py export vcf $(inputs.output_basename).call.cns -o $(inputs.output_basename).vcf

      cnvkit.py metrics $(inputs.input_sample.nameroot).cnr -s $(inputs.input_sample.nameroot).cns
      -o $(inputs.output_basename).metrics.txt

      cnvkit.py gainloss --drop-low-coverage $(inputs.input_sample.nameroot).cnr -o $(inputs.output_basename).gainloss.txt

      mv $(inputs.input_sample.nameroot).cnr $(inputs.output_basename).cnr

      mv $(inputs.input_sample.nameroot)-diagram.pdf $(inputs.output_basename).diagram.pdf
      
      mv $(inputs.input_sample.nameroot)-scatter.pdf $(inputs.output_basename).scatter.pdf


inputs:
  input_sample: {type: File, doc: "tumor bam file", secondaryFiles: [^.bai]}
  input_control: {type: ['null', File], doc: "normal bam file", secondaryFiles: [^.bai]}
  reference: {type: File, doc: "fasta file", secondaryFiles: [.fai]}
  cnv_kit_cnn: {type: ['null', File], doc: "If running using an existing .cnn, supply here"}
  b_allele_vcf: {type: ['null', File], doc: "b allele germline vcf, if available"}
  capture_regions: {type: ['null', File], doc: "target regions for WES"}
  annotation_file: {type: File, doc: "refFlat.txt file"}
  output_basename: string
  wgs_mode: {type: ['null', string], doc: "for WGS mode, input Y. leave blank for hybrid mode"}
  threads:
    type: ['null', int]
    default: 16

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
      glob: '*.cnn'
    doc: "Output if starting from cnn scratch.  Should not appear if an existing .cnn was given as input."
