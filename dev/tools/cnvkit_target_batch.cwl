cwlVersion: v1.0
class: CommandLineTool
id: cnvkit-rs-target-batch
doc: "This is a research version of the prod tool, fir targeted only"
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
      $(inputs.input_sample.path)
      --targets $(inputs.capture_regions.path)
      ${
        if (inputs.avg_target_size){
          return "--target-avg-size " + inputs.avg_target_size;
        }
        else{
          return "";
        }
      }
      ${
        if (inputs.drop_low_coverage){
          return "--drop-low-coverage"
        }
        else{
          return "";
        }
      }
      ${
        if(inputs.antitargets){
          return "--antitargets " +  inputs.antitargets.path
        }
        else{
          return "";
        }
      }
      ${
        if (inputs.cnv_kit_cnn == null){
          var arg = "--output-reference " + inputs.output_basename + "_cnvkit_reference.cnn --fasta " + inputs.reference.path + " --annotate " + inputs.annotation_file.path;
          if (inputs.input_control != null) {
              arg += " --normal " + inputs.input_control.path;
          }
        }
        else{
          var arg = "--reference " + inputs.cnv_kit_cnn.path;
          var msex = ['m','y','male','Male']
          if (msex.indexOf(inputs.sex) >= 0){
            arg += " --male-reference";
          }
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
      ${
        var arg = "--sample-sex " + inputs.sex;
        var msex = ['m','y','male','Male']
        if (msex.indexOf(inputs.sex) >= 0){
          arg += " --male-reference";
        }
        return arg;
      }
      -o $(inputs.output_basename).call.cns
            
      ln -s $(inputs.output_basename).call.cns $(inputs.tumor_sample_name).cns

      cnvkit.py export seg $(inputs.tumor_sample_name).cns -o $(inputs.output_basename).call.seg

      rm $(inputs.tumor_sample_name).cns

      cnvkit.py metrics $(inputs.input_sample.nameroot).cnr -s $(inputs.input_sample.nameroot).cns
      -o $(inputs.output_basename).metrics.txt

      cnvkit.py gainloss $(inputs.input_sample.nameroot).cnr -o $(inputs.output_basename).gainloss.txt

      mv $(inputs.input_sample.nameroot).cnr $(inputs.output_basename).cnr

      mv $(inputs.input_sample.nameroot)-diagram.pdf $(inputs.output_basename).diagram.pdf
      
      mv $(inputs.input_sample.nameroot)-scatter.pdf $(inputs.output_basename).scatter.pdf


inputs:
  input_sample: {type: File, doc: "tumor bam file", secondaryFiles: [^.bai]}
  input_control: {type: ['null', File], doc: "normal bam file", secondaryFiles: [^.bai]}
  reference: {type: ['null', File], doc: "fasta file, needed if cnv kit cnn not already built", secondaryFiles: [.fai]}
  cnvkit_cnn: {type: ['null', File], doc: "If running using an existing .cnn, supply here"}
  b_allele_vcf: {type: ['null', File], doc: "b allele germline vcf, if available"}
  capture_regions: {type: ['null', File], doc: "target regions for WES"}
  annotation_file: {type: ['null', File], doc: "refFlat.txt file,  needed if cnv kit cnn not already built"}
  output_basename: string
  tumor_sample_name: string
  threads:
    type: ['null', int]
    default: 16
  sex:
    type: string
    doc: "Set sample sex.  CNVkit isn't always great at guessing it"
  avg_target_size: {type: int?, doc: "CNVkit recommends 267, smaller can be used but urged to also have anti-target file"}
  drop_low_coverage: {type: boolean, doc: "Drop very-low-coverage bins before segmentation to avoid false-positive deletions in poor-quality tumor samples", default: false}
  antitargets: {type: File?, doc: "A bed file with anti-targets"}


outputs:
  output_cnr: 
    type: File
    outputBinding:
      glob: '*.cnr'
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