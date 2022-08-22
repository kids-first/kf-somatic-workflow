cwlVersion: v1.0
class: CommandLineTool
id: kf-mskcc-vcf2maf-1.6.21
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 4000
    coresMin: 2
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/kf_vcf2maf:v1.0.3'
  - class: InitialWorkDirRequirement
    listing:
      - entryname: vcf2maf.pl
        entry:
          $include: ../scripts/kf_mskcc_vcf2maf_1.6.21.pl

baseCommand: [gunzip, -c]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.input_vcf.path) > input_file.vcf
      && perl vcf2maf.pl
      --input-vcf input_file.vcf
      --output-maf $(inputs.output_basename).$(inputs.tool_name).vep.maf
      ${
        if(inputs.maf_center){
          return "--maf-center \"" + inputs.maf_center + "\""
        }
        else{
          return "";
        }
      }

inputs:
  reference: { type: 'File',  secondaryFiles: [.fai], doc: "Fasta genome assembly with index",
    inputBinding: {position: 2, prefix: "--ref-fasta"} }
  input_vcf: { type: 'File', secondaryFiles: [.tbi], doc: "VEP annotated vcf file." }
  output_basename: string
  tumor_id: { type: string, inputBinding: {position: 3, prefix: "--tumor-id"} }
  normal_id: { type: string, inputBinding: {position: 4, prefix: "--normal-id"} }
  tool_name: string
  ref_build: { type: 'string?', doc: "Genome ref build used, should line up with cache.", default: "GRCh38",
    inputBinding: {position: 5, prefix: "--ncbi-build"} }
  retain_info: { type: 'string?', doc: "csv string with INFO fields that you want to keep, i.e. for consensus `MQ,MQ0,CAL,Hotspot`",
    inputBinding: {position: 6, prefix: "--retain-info"} }
  retain_fmt: { type: 'string?', doc: "csv string with FORMAT fields that you want to keep",
    inputBinding: {position: 6, prefix: "--retain-fmt"} }
  retain_ann: { type: 'string?', doc: "csv string of annotations (within the VEP CSQ/ANN) to retain as extra columns in MAF",
    inputBinding: {position: 6, prefix: "--retain-ann"} }
  custom_enst: { type: 'File?', doc: "Use a file with ens tx IDs for each gene to override VEP PICK",
    inputBinding: {position: 7, prefix: "--custom-enst"} }
  maf_center: {type: 'string?', doc: "Sequencing center of variant called", default: "."}

outputs:
  output_maf:
    type: File
    outputBinding:
      glob: '*.maf'
