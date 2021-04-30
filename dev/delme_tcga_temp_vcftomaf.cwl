cwlVersion: v1.0
class: Workflow
id: tcga_vcf2maf_patch
doc: "This is a temp wf to convert vcfs to maf and rename"
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  reference: {type: File,  secondaryFiles: [.fai], doc: "Fasta genome assembly with index"}
  strelka2_public_vcf: {type: File, secondaryFiles: [.tbi], doc: "VEP annotated vcf file."}
  mutect2_public_vcf: {type: File, secondaryFiles: [.tbi], doc: "VEP annotated vcf file."}
  lancet_public_vcf: {type: File, secondaryFiles: [.tbi], doc: "VEP annotated vcf file."}
  strelka2_protected_vcf: {type: File, secondaryFiles: [.tbi], doc: "VEP annotated vcf file."}
  mutect2_protected_vcf: {type: File, secondaryFiles: [.tbi], doc: "VEP annotated vcf file."}
  lancet_protected_vcf: {type: File, secondaryFiles: [.tbi], doc: "VEP annotated vcf file."}
  output_basename: string
  tumor_id: string
  normal_id: string

outputs:
  protected_mafs: {type: 'File[]', outputSource: rename_outputs_protected/renamed_files}
  public_mafs: {type: 'File[]', outputSource: rename_outputs_public/renamed_files}

steps:
  strelka2_vcf2maf_protected:
    run: ../tools/kf_mskcc_vcf2maf.cwl
    in:
      reference: reference
      input_vcf: strelka2_protected_vcf
      output_basename: output_basename
      tumor_id: tumor_id
      normal_id: normal_id
      tool_name:
        valueFrom: ${return "strelka2_somatic"}
      retain_info:
        valueFrom: ${return "MQ,MQ0,QSI,HotSpotAllele"}
    out: [output_maf]

  mutect2_vcf2maf_protected:
    run: ../tools/kf_mskcc_vcf2maf.cwl
    in:
      reference: reference
      input_vcf: mutect2_protected_vcf
      output_basename: output_basename
      tumor_id: tumor_id
      normal_id: normal_id
      tool_name:
        valueFrom: ${return "mutect2_somatic"}
      retain_info:
        valueFrom: ${return "MBQ,TLOD,HotSpotAllele"}
    out: [output_maf]

  lancet_vcf2maf_protected:
    run: ../tools/kf_mskcc_vcf2maf.cwl
    in:
      reference: reference
      input_vcf: lancet_protected_vcf
      output_basename: output_basename
      tumor_id: tumor_id
      normal_id: normal_id
      tool_name:
        valueFrom: ${return "lancet_somatic"}
      retain_info:
        valueFrom: ${return "MS,FETS,HotSpotAllele"}
    out: [output_maf]

  strelka2_vcf2maf_public:
    run: ../tools/kf_mskcc_vcf2maf.cwl
    in:
      reference: reference
      input_vcf: strelka2_public_vcf
      output_basename: output_basename
      tumor_id: tumor_id
      normal_id: normal_id
      tool_name:
        valueFrom: ${return "strelka2_somatic"}
      retain_info:
        valueFrom: ${return "MQ,MQ0,QSI,HotSpotAllele"}
    out: [output_maf]

  mutect2_vcf2maf_public:
    run: ../tools/kf_mskcc_vcf2maf.cwl
    in:
      reference: reference
      input_vcf: mutect2_public_vcf
      output_basename: output_basename
      tumor_id: tumor_id
      normal_id: normal_id
      tool_name:
        valueFrom: ${return "mutect2_somatic"}
      retain_info:
        valueFrom: ${return "MBQ,TLOD,HotSpotAllele"}
    out: [output_maf]

  lancet_vcf2maf_public:
    run: ../tools/kf_mskcc_vcf2maf.cwl
    in:
      reference: reference
      input_vcf: lancet_public_vcf
      output_basename: output_basename
      tumor_id: tumor_id
      normal_id: normal_id
      tool_name:
        valueFrom: ${return "lancet_somatic"}
      retain_info:
        valueFrom: ${return "MS,FETS,HotSpotAllele"}
    out: [output_maf]

  rename_outputs_protected:
    run: ../tools/generic_rename_outputs.cwl
    in:
      input_files:
        source: [strelka2_vcf2maf_protected/output_maf, mutect2_vcf2maf_protected/output_maf, lancet_vcf2maf_protected/output_maf]
        valueFrom: "${return [ self[0],self[1],self[2] ]}"
      rename_to:
        source: output_basename
        valueFrom: "${var str_maf=self + '.strelka2_somatic.norm.annot.protected.maf'; \
        var mut_maf=self + '.mutect2_somatic.norm.annot.protected.maf'; \
        var lan_maf=self + '.lancet_somatic.norm.annot.protected.maf'; \
        return [str_maf, mut_maf, lan_maf];}"
    out: [renamed_files]

  rename_outputs_public:
    run: ../tools/generic_rename_outputs.cwl
    in:
      input_files:
        source: [strelka2_vcf2maf_public/output_maf, mutect2_vcf2maf_public/output_maf, lancet_vcf2maf_public/output_maf]
        valueFrom: "${return [ self[0],self[1],self[2] ]}"
      rename_to:
        source: output_basename
        valueFrom: "${var str_maf=self + '.strelka2_somatic.norm.annot.public.maf'; \
        var mut_maf=self + '.mutect2_somatic.norm.annot.public.maf'; \
        var lan_maf=self + '.lancet_somatic.norm.annot.public.maf'; \
        return [str_maf, mut_maf, lan_maf];}"
    out: [renamed_files]
