cwlVersion: v1.0
class: Workflow
id: kfdrc_consensus_calling
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: ['.fai', '^.dict']}
  strelka2_vcf: {type: File, secondaryFiles: ['.tbi']}
  mutect2_vcf: {type: File, secondaryFiles: ['.tbi']}
  lancet_vcf: {type: File, secondaryFiles: ['.tbi']}
  vardict_vcf: {type: File, secondaryFiles: ['.tbi']}
  cram: {type: File, secondaryFiles: ['.crai']}
  input_tumor_name: string
  input_normal_name: string
  output_basename: string
  ncallers: {type: int?, doc: "Optional number of callers required for consensus [2]", default: 2}
  annotation_vcf: {type: File, secondaryFiles: ['.tbi'], doc: "VCF of annotions to add to consensus variants, e.g. gnomAD allele frequency"}
  vep_cache: {type: File, doc: "tar gzipped cache from ensembl/local converted cache"}
  annot_columns: {type: string?, default: 'INFO/AF', doc: "column from annotation_vcf to add to consensus VCF; defaults to 'INFO/AF'"}
  filter_names: {type: 'string[]?', default: [ "NORM_DP_LOW", "GNOMAD_AF_HIGH" ], doc: "Names of filters to be added to consensus VCF;\
     \ default values set"}
  depth_lowerbound: {type: int?, default: 7, doc: "Normal-sample read depth at which to apply depth filter; default set"}
  frequency_upperbound: {type: float?, default: 0.001, doc: "Population allele frequency above which to apply frequency filter; default set"}
  bcftools_public_filter: {type: string?, doc: 'Will hard filter final result to create a public version, e.g. FILTER="PASS"|INFO/Hotspot=1; default set', default: 'FILTER="PASS"|INFO/Hotspot=1'}
  retain_info: {type: string?, doc: "csv string with INFO fields that you want to keep; default values set", default: 'MQ,MQ0,CAL,HotSpotAllele'}
  retain_fmt: {type: string?, doc: "csv string with FORMAT fields that you want to keep"}
  use_kf_fields: {type: boolean?, doc: "Flag to drop fields normally not used in KF, or keep cBio defaults", default: true}

outputs:
  annotated_protected_outputs: {type: 'File[]', outputSource: rename_protected/renamed_files}
  annotated_public_outputs: {type: 'File[]', outputSource: rename_public/renamed_files}

steps:
  prep_mnp_variants:
    run: ../tools/prep_mnp_variants.cwl
    in:
      strelka2_vcf: strelka2_vcf
      other_vcfs: [mutect2_vcf, lancet_vcf, vardict_vcf]
      output_basename: output_basename
    out: [output_vcfs]

  consensus_merge:
    run: ../tools/consensus_merge.cwl
    in:
      strelka2_vcf:
        source: prep_mnp_variants/output_vcfs 
        valueFrom: '$(self[0])'
      mutect2_vcf: mutect2_vcf
      lancet_vcf: lancet_vcf
      vardict_vcf: vardict_vcf
      cram: cram
      ncallers: ncallers
      reference: indexed_reference_fasta
      output_basename: output_basename
    out: [output]

  vep_annot_consensus:
    run: ../tools/vep_somatic_annotate_r93.cwl
    in:
      input_vcf: consensus_merge/output
      output_basename: output_basename
      tool_name:
        valueFrom: ${return "consensus_somatic"}
      reference: indexed_reference_fasta
      cache: vep_cache
    out: [output_vcf]

  bcftools_annotate:
    run: ../tools/bcftools_annotate.cwl
    in:
      input_vcf: vep_annot_consensus/output_vcf
      annotation_vcf: annotation_vcf
      output_basename: output_basename
      tool_name:
        valueFrom: ${return "bcft_annot"}
      columns: annot_columns
    out: [bcftools_annotated_vcf]

  variant_filter:
    run: ../tools/gatk_variant_filter.cwl
    in:
      input_vcf: bcftools_annotate/bcftools_annotated_vcf
      output_basename: output_basename
      reference: indexed_reference_fasta
      tool_name:
        valueFrom: ${return "filter"}
      filter_name: filter_names
      filter_expression:
        source: [input_normal_name, depth_lowerbound, frequency_upperbound]
        valueFrom: ${return [ "vc.getGenotype(" + self[0] + ").getDP() <= " + self[1], "AF > " + self[2] ]}
    out: [gatk_soft_filtered_vcf]

  vcf2maf:
    run: ../tools/kf_mskcc_vcf2maf.cwl
    in:
      input_vcf: variant_filter/gatk_soft_filtered_vcf
      normal_id: input_normal_name
      tumor_id: input_tumor_name
      output_basename: output_basename
      reference: indexed_reference_fasta
      retain_info: retain_info
      tool_name:
        valueFrom: ${return "vcf2maf"}
    out: [output_maf]

  hard_filter_vcf:
    run: ../tools/bcftools_filter_vcf.cwl
    in:
      input_vcf: variant_filter/gatk_soft_filtered_vcf
      include_expression: bcftools_public_filter
      output_basename: output_basename
    out:
      [filtered_vcf]

  kfdrc_vcf2maf_public:
    run: ../tools/kf_mskcc_vcf2maf.cwl
    in:
      reference: indexed_reference_fasta
      input_vcf: hard_filter_vcf/filtered_vcf
      output_basename: output_basename
      tumor_id: input_tumor_name
      normal_id: input_normal_name
      tool_name:
        valueFrom: ${return "public"}
      retain_info: retain_info
      retain_fmt: retain_fmt
      use_kf_fields: use_kf_fields
    out: [output_maf]

  rename_protected:
    run: ../tools/generic_rename_outputs.cwl
    in:
      input_files:
        source: [variant_filter/gatk_soft_filtered_vcf, vcf2maf/output_maf]
        valueFrom: "${return [self[0],self[0].secondaryFiles[0],self[1]]}"
      rename_to:
        source: output_basename
        valueFrom: "${var pro_vcf=self + '.protected.vcf.gz'; \
        var pro_tbi=self + '.protected.vcf.gz.tbi'; \
        var pro_maf=self + '.protected.maf'; \
        return [pro_vcf, pro_tbi, pro_maf];}"
    out: [renamed_files]

  rename_public:
    run: ../tools/generic_rename_outputs.cwl
    in:
      input_files:
        source: [hard_filter_vcf/filtered_vcf, kfdrc_vcf2maf_public/output_maf]
        valueFrom: "${return [self[0],self[0].secondaryFiles[0],self[1]]}"
      rename_to:
        source: output_basename
        valueFrom: "${var pub_vcf=self + '.public.vcf.gz'; \
        var pub_tbi=self + '.public.vcf.gz.tbi'; \
        var pub_maf=self + '.public.maf'; \
        return [pub_vcf, pub_tbi, pub_maf];}"
    out: [renamed_files]
