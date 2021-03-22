cwlVersion: v1.0
class: Workflow
id: kfdrc_consensus_calling
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai]}
  strelka2_vcf: {type: File, secondaryFiles: ['.tbi']}
  mutect2_vcf: {type: File, secondaryFiles: ['.tbi']}
  lancet_vcf: {type: File, secondaryFiles: ['.tbi']}
  vardict_vcf: {type: File, secondaryFiles: ['.tbi']}
  cram: {type: File, secondaryFiles: ['.crai']}
  input_tumor_name: string
  input_normal_name: string
  output_basename: string
  annotation_vcf: {type: File, secondaryFiles: ['.tbi'], doc: "VCF of annotions to add to consensus variants, e.g. gnomAD allele frequency"}
  vep_cache: {type: File, doc: "tar gzipped cache from ensembl/local converted cache"}
  strip_info: {type: ['null', string], doc: "If given, remove previous annotation information based on INFO file, i.e. to strip VEP info, use INFO/CSQ or INFO/ANN depending on how VEP was run for instance"}

outputs:
  vep_consensus_vcf: {type: File, outputSource: variant_filter/gatk_soft_filtered_vcf, secondaryFiles: ['.tbi']}
  vep_consensus_maf: {type: File, outputSource: vcf2maf/output_maf}

steps:
  prep_mnp_variants:
    run: ../tools/prep_mnp_variants.cwl
    in:
      strelka2_vcf: strelka2_vcf
      other_vcfs: [mutect2_vcf, lancet_vcf, vardict_vcf]
      output_basename: output_basename
    out: [output_vcfs]

  consensus_merge:
    hints:
    - class: 'sbg:AWSInstanceType'
      value: c5.4xlarge
    run: ../tools/consensus_merge.cwl
    in:
      strelka2_vcf:
        source: prep_mnp_variants/output_vcfs 
        valueFrom: '$(self[0])'
      mutect2_vcf: mutect2_vcf
      lancet_vcf: lancet_vcf
      vardict_vcf: vardict_vcf
      cram: cram
      reference: indexed_reference_fasta
      output_basename: output_basename
    out: [output]

  vep_annot_consensus:
    run: ../tools/vep_somatic_annotate_r93.cwl
    in:
      input_vcf: consensus_merge/output
      output_basename: output_basename
      tumor_id: input_tumor_name
      normal_id: input_normal_name
      tool_name:
        valueFrom: ${return "consensus_somatic"}
      reference: indexed_reference_fasta
      cache: vep_cache
    out: [output_vcf, output_tbi]

  bcftools_annotate:
    run: ../tools/bcftools_annotate.cwl
    in:
      input_vcf: vep_annot_consensus/output_vcf
      annotation_vcf: annotation_vcf
      output_basename: output_basename
      tool_name:
        valueFrom: ${return "bcft_annot"}
      columns:
        valueFrom: ${return "INFO/AF"}
    out: [bcftools_annotated_vcf]

  variant_filter:
    run: ../tools/gatk_variant_filter.cwl
    in:
      input_vcf: bcftools_annotate/bcftools_annotated_vcf
      output_basename: output_basename
      reference: indexed_reference_fasta
      tool_name:
        valueFrom: ${return "filter"}
      filter_name:
        valueFrom: ${return [ "NORM_DP_LOW", "GNOMAD_AF_HIGH" ]}
      filter_expression:
        valueFrom: ${return [ "vc.getGenotype(<normal sample name>).getDP() <= 7", "AF > 0.001" ]}
    out: [gatk_soft_filtered_vcf]

  vcf2maf:
    run: ../tools/kf_mskcc_vcf2maf.cwl
    in:
      input_vcf: variant_filter/gatk_soft_filtered_vcf
      normal_id: input_normal_name
      tumor_id: input_tumor_name
      output_basename: output_basename
      reference: indexed_reference_fasta
      retain_info:
        valueFrom: ${return [ "MQ" ,"MQ0", "CAL", "HotSpotAllele" ]}
      tool_name:
        valueFrom: ${return "vcf2maf"}
    out: [output_maf]

$namespaces:
  sbg: https://sevenbridges.com
