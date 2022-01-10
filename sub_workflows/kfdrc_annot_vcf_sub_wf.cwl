cwlVersion: v1.0
class: Workflow
id: kfdrc_annot_sub_wf
doc: "This subworkflow normalizes, annotates, and adds soft filters to vcf inputs"
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  indexed_reference_fasta: {type: 'File', secondaryFiles: [.fai, ^.dict]}
  input_vcf: {type: 'File', secondaryFiles: ['.tbi'], doc: "Input vcf to annotate and soft filter"}
  input_tumor_name: string
  input_normal_name: string
  add_common_fields: {type: 'boolean', doc: "Set to true if input is a strelka2 vcf that hasn't had common fields added", default: false}
  bcftools_annot_columns: {type: 'string', doc: "csv string of columns from annotation to port into the input vcf, i.e INFO/AF"}
  bcftools_annot_vcf: {type: 'File', secondaryFiles: ['.tbi'], doc: "bgzipped annotation vcf file"}
  bcftools_public_filter: {type: 'string?', doc: "Will hard filter final result to create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  gatk_filter_name: {type: 'string[]', doc: "Array of names for each filter tag to add"}
  gatk_filter_expression: {type: 'string[]', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration for clues"}
  vep_cache: {type: 'File', doc: "tar gzipped cache from ensembl/local converted cache"}
  vep_ref_build: {type: ['null', string], doc: "Genome ref build used, should line up with cache.", default: "GRCh38" }
  disable_hotspot_annotation: { type: 'boolean?', doc: "Disable Hotspot Annotation and skip this task." }
  genomic_hotspots: { type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots" }
  protein_snv_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots" }
  protein_indel_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspots" }
  output_basename: string
  tool_name: string
  retain_info: {type: 'string?', doc: "csv string with INFO fields that you want to keep, i.e. for consensus `MQ,MQ0,CAL,Hotspot`"}
  retain_fmt: {type: 'string?', doc: "csv string with FORMAT fields that you want to keep"}
  maf_center: {type: 'string?', doc: "Sequencing center of variant called", default: "."}

outputs:
  annotated_protected_vcf: {type: 'File', outputSource: hotspots_annotation/hotspots_vcf}
  annotated_protected_maf: {type: 'File', outputSource: kfdrc_vcf2maf_protected/output_maf}
  annotated_public_vcf: {type: 'File', outputSource: hard_filter_vcf/filtered_vcf}
  annotated_public_maf: {type: 'File', outputSource: kfdrc_vcf2maf_public/output_maf}

steps:
  normalize_vcf:
    run: ../tools/normalize_vcf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf: input_vcf
      output_basename: output_basename
      tool_name: tool_name
    out: [normalized_vcf]

  bcftools_strip_info:
    run: ../tools/bcftools_strip_ann.cwl
    in:
      input_vcf: normalize_vcf/normalized_vcf
      output_basename: output_basename
      tool_name: tool_name
      strip_info: bcftools_annot_columns
    out: [stripped_vcf]

  add_standard_fields:
    run: ../tools/add_strelka2_fields.cwl
    in:
      strelka2_vcf: bcftools_strip_info/stripped_vcf
      run_tool_flag: add_common_fields
      tumor_name: input_tumor_name
      normal_name: input_normal_name
      output_basename: output_basename
    out: [output]

  vep_annotate_vcf:
    run: ../tools/vep_somatic_annotate_r93.cwl
    in:
      reference: indexed_reference_fasta
      input_vcf: add_standard_fields/output
      output_basename: output_basename
      tool_name: tool_name
      cache: vep_cache
      ref_build: vep_ref_build
    out: [output_vcf]

  bcftools_gnomad_annotate:
    run: ../tools/bcftools_annotate.cwl
    in:
      input_vcf: vep_annotate_vcf/output_vcf
      annotation_vcf: bcftools_annot_vcf
      columns: bcftools_annot_columns
      output_basename: output_basename
      tool_name: tool_name
    out: [bcftools_annotated_vcf]

  gatk_add_soft_filter:
    run: ../tools/gatk_variant_filter.cwl
    in:
      input_vcf: bcftools_gnomad_annotate/bcftools_annotated_vcf
      reference: indexed_reference_fasta
      filter_name: gatk_filter_name
      filter_expression: gatk_filter_expression
      output_basename: output_basename
      tool_name: tool_name
    out: [gatk_soft_filtered_vcf]

  hotspots_annotation:
    run: ../tools/hotspots_annotation.cwl
    in:
      input_vcf: gatk_add_soft_filter/gatk_soft_filtered_vcf
      disable_hotspot_annotation: disable_hotspot_annotation
      genomic_hotspots: genomic_hotspots
      protein_snvs: protein_snv_hotspots
      protein_indels: protein_indel_hotspots
      output_basename: output_basename
    out: [hotspots_vcf]

  kfdrc_vcf2maf_protected:
    run: ../tools/kf_mskcc_vcf2maf.cwl
    in:
      reference: indexed_reference_fasta
      input_vcf: hotspots_annotation/hotspots_vcf
      output_basename: output_basename
      tumor_id: input_tumor_name
      normal_id: input_normal_name
      tool_name: tool_name
      retain_info: retain_info
      retain_fmt: retain_fmt
      maf_center: maf_center
    out: [output_maf]

  hard_filter_vcf:
    run: ../tools/bcftools_filter_vcf.cwl
    in:
      input_vcf: hotspots_annotation/hotspots_vcf
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
      tool_name: tool_name
      retain_info: retain_info
      retain_fmt: retain_fmt
      maf_center: maf_center
    out: [output_maf]
