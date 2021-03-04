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
    dockerPull: 'migbro/vcf2maf:v1.0.1'
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
      --custom-enst /vcf2maf/data/isoform_overrides_uniprot
      ${
        if(inputs.use_kf_fields){
          return "--use-kf-fields";
        }
        else{
          return "";
        }
      }
      --ncbi-build $(inputs.ref_build)
      ${
        if(inputs.retain_info){
          return "--retain-info " + inputs.retain_info;
        }
        else{
          return "";
        }
      }
      --ref-fasta $(inputs.reference.path)

inputs:
  reference: {type: File,  secondaryFiles: [.fai], doc: "Fasta genome assembly with index"}
  input_vcf: {type: File, secondaryFiles: [.tbi], doc: "VEP annotated vcf file."}
  output_basename: string
  tumor_id: string
  normal_id: string
  tool_name: string
  ref_build: {type: string?, doc: "Genome ref build used, should line up with cache.", default: "GRCh38"}
  retain_info: {type: string?, doc: "csv string with INFO fields that you want to keep, i.e. for consensus `MQ,MQ0,CAL,Hotspot`"}
  use_kf_fields: {type: boolean?, doc: "Flag to drop fields normally not used in KF, or keep cBio defaults", default: true}

outputs:
  output_maf:
    type: File
    outputBinding:
      glob: '*.maf'