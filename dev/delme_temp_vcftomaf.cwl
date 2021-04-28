cwlVersion: v1.0
class: Workflow
id: patch_vcf2maf
doc: "This is a temp wf to convert vcfs to maf and rename"
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  reference: {type: File,  secondaryFiles: [.fai], doc: "Fasta genome assembly with index"}
  public_vcf: {type: File, secondaryFiles: [.tbi], doc: "VEP annotated vcf file."}
  protected_vcf: {type: File, secondaryFiles: [.tbi], doc: "VEP annotated vcf file."}
  output_basename: string
  tumor_id: string
  normal_id: string
  tool_name: string
  ref_build: {type: string?, doc: "Genome ref build used, should line up with cache.", default: "GRCh38"}
  retain_info: {type: string?, doc: "csv string with INFO fields that you want to keep, i.e. for consensus `MQ,MQ0,CAL,Hotspot`"}
  retain_fmt: {type: string?, doc: "csv string with FORMAT fields that you want to keep"}
  maf_center: {type: string?, doc: "Sequencing center of variant called", default: "."}

outputs:
  patched_outputs: {type: 'File[]', outputSource: rename_outputs/renamed_files}

steps:
  kfdrc_vcf2maf_protected:
    run: ../tools/kf_mskcc_vcf2maf.cwl
    in:
      reference: reference
      input_vcf: protected_vcf
      output_basename: output_basename
      tumor_id: tumor_id
      normal_id: normal_id
      tool_name: tool_name
      retain_info: retain_info
      retain_fmt: retain_fmt
      ref_build: ref_build
      maf_center: maf_center
    out: [output_maf]

  kfdrc_vcf2maf_public:
    run: ../tools/kf_mskcc_vcf2maf.cwl
    in:
      reference: reference
      input_vcf: public_vcf
      output_basename: output_basename
      tumor_id: tumor_id
      normal_id: normal_id
      tool_name: tool_name
      retain_info: retain_info
      retain_fmt: retain_fmt
      ref_build: ref_build
      maf_center: maf_center
    out: [output_maf]

  rename_outputs:
    run: ../tools/generic_rename_outputs.cwl
    in:
      input_files:
        source: [kfdrc_vcf2maf_protected/output_maf, kfdrc_vcf2maf_public/output_maf]
        valueFrom: "${return [ self[0],self[1] ]}"
      rename_to:
        source: [output_basename, tool_name]
        valueFrom: "${var pro_maf=self[0] + '.' + self[1] + '.norm.annot.protected.maf'; \
        var pub_maf=self[0] + '.' + self[1] + '.norm.annot.public.maf'; \
        return [pro_maf, pub_maf];}"
    out: [renamed_files]

