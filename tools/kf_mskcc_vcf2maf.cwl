cwlVersion: v1.0
class: CommandLineTool
id: kf-mskcc-vcf2maf
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 4000
    coresMin: 2
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/kf_vcf2maf:v1.0.3'
baseCommand: [gunzip, -c]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.input_vcf.path) > input_file.vcf
      && perl /vcf2maf/vcf2maf.pl
      --input-vcf input_file.vcf
      --output-maf $(inputs.output_basename).$(inputs.tool_name).vep.maf
      --tumor-id $(inputs.tumor_id)
      --normal-id $(inputs.normal_id)
      --ncbi-build $(inputs.ref_build)
      --ref-fasta $(inputs.reference.path)
      ${
        if(inputs.maf_center){
          return "--maf-center \"" + inputs.maf_center + "\""
        }
        else{
          return "";
        }
      }
      ${
        if(inputs.retain_info){
          return "--retain-info " + inputs.retain_info;
        }
        else{
          return "";
        }
      }
      ${
        if(inputs.retain_fmt){
          return "--retain-fmt " + inputs.retain_fmt;
        }
        else{
          return "";
        }
      }
      ${
        if(inputs.custom_enst){
          return "--custom-enst " + inputs.custom_enst.path;
        }
        else{
          return "";
        }
      }

inputs:
  reference: {type: File,  secondaryFiles: [.fai], doc: "Fasta genome assembly with index"}
  input_vcf: {type: File, secondaryFiles: [.tbi], doc: "VEP annotated vcf file."}
  output_basename: string
  tumor_id: string
  normal_id: string
  tool_name: string
  ref_build: {type: string?, doc: "Genome ref build used, should line up with cache.", default: "GRCh38"}
  retain_info: {type: string?, doc: "csv string with INFO fields that you want to keep, i.e. for consensus `MQ,MQ0,CAL,Hotspot`"}
  retain_fmt: {type: string?, doc: "csv string with FORMAT fields that you want to keep"}
  custom_enst: {type: File?, doc: "Use a file with ens tx IDs for each gene to override VEP PICK"}
  maf_center: {type: string?, doc: "Sequencing center of variant called", default: "."}

outputs:
  output_maf:
    type: File
    outputBinding:
      glob: '*.maf'
