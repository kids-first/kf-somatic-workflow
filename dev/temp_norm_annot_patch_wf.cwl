cwlVersion: v1.0
class: Workflow
id: kfdrc_norm_annot_wf
doc: "This is a temp wf to normalize and annotate a PASS vcf"
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict]}
  input_vcf: {type: File, secondaryFiles: ['.tbi'], doc: "Input vcf to normalize, annotate, and soft filter"}
  input_tumor_name: string
  input_normal_name: string
  add_common_fields: {type: boolean, doc: "Set to true if input is a strelka2 vcf that hasn't had common fields added", default: false}
  bcftools_annot_columns: {type: string, doc: "csv string of columns from annotation to port into the input vcf, i.e INFO/AF", default: "INFO/AF"}
  bcftools_annot_vcf: {type: File, secondaryFiles: ['.tbi'], doc: "bgzipped annotation vcf file"}
  gatk_filter_name: {type: 'string[]', doc: "Array of names for each filter tag to add, recommend: [\"NORM_DP_LOW\", \"GNOMAD_AF_HIGH\"]"}
  gatk_filter_expression: {type: 'string[]', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration, recommend: \"vc.getGenotype('\" + inputs.input_normal_name + \"').getDP() <= 7\"), \"AF > 0.001\"]"}
  disable_hotspot_annotation: { type: 'boolean?', doc: "Disable Hotspot Annotation and skip this task.", default: false }
  vep_cache: {type: File, doc: "tar gzipped cache from ensembl/local converted cache"}
  vep_ref_build: {type: ['null', string], doc: "Genome ref build used, should line up with cache.", default: "GRCh38" }
  genomic_hotspots: { type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots" }
  protein_snv_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots" }
  protein_indel_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspots" }
  output_basename: string
  tool_name: string
  retain_info: {type: string?, doc: "csv string with INFO fields that you want to keep, i.e. for consensus `MQ,MQ0,CAL,HotSpotAllele`"}
  retain_fmt: {type: string?, doc: "csv string with FORMAT fields that you want to keep"}
  use_kf_fields: {type: boolean?, doc: "Flag to drop fields normally not used in KF, or keep cBio defaults", default: true}

outputs:
  annotated_output: {type: 'File[]', outputSource: rename_outputs/renamed_files}

steps:
  normalize_vcf:
    run: ../tools/normalize_vcf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf: input_vcf
      output_basename: output_basename
      tool_name: tool_name
    out: [normalized_vcf]

  annotate:
    run: ../sub_workflows/kfdrc_annot_vcf_sub_wf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf: normalize_vcf/normalized_vcf
      input_tumor_name: input_tumor_name
      input_normal_name: input_normal_name
      add_common_fields: add_common_fields
      retain_info: retain_info
      use_kf_fields: use_kf_fields
      bcftools_annot_columns: bcftools_annot_columns
      bcftools_annot_vcf: bcftools_annot_vcf
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      vep_cache: vep_cache
      vep_ref_build: vep_ref_build
      disable_hotspot_annotation: disable_hotspot_annotation
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      output_basename: output_basename
      tool_name: tool_name
    out: [annotated_vcf, annotated_maf]

  rename_outputs:
    run: ../tools/generic_rename_outputs.cwl
    in:
      input_files:
        source: [annotate/annotated_vcf, annotate/annotated_maf]
        valueFrom: "${return [self[0], self[0].secondaryFiles[0], self[1]]}"
      rename_to:
        source: [output_basename, tool_name]
        valueFrom: "${var vcf=self[0] + '.' + self[1] + '.norm.annot.vcf.gz'; var tbi = self[0] + '.' + self[1] + '.norm.annot.vcf.gz.tbi'; var maf = self[0] + '.' + self[1] + '.norm.annot.maf'; return [vcf, tbi, maf];}"
    out: [renamed_files]


