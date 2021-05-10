cwlVersion: v1.0
class: Workflow
id: kfdrc_production_strelka2_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
inputs:
  # Required
  reference_fasta: { type: File, sbg:suggestedValue: {class: File, path: 60639014357c3a53540ca7a3, name: Homo_sapiens_assembly38.fasta} }
  reference_fai: { type: 'File?', sbg:suggestedValue: {class: File, path: 60639016357c3a53540ca7af, name: Homo_sapiens_assembly38.fasta.fai} }
  reference_dict: { type: 'File?', sbg:suggestedValue: {class: File, path: 60639019357c3a53540ca7e7, name: Homo_sapiens_assembly38.dict} }
  input_tumor_aligned:
    type: File
    secondaryFiles: |
      ${
        var dpath = self.location.replace(self.basename, "")
        if(self.nameext == '.bam'){
          return {"location": dpath+self.nameroot+".bai", "class": "File"}
        }
        else{
          return {"location": dpath+self.basename+".crai", "class": "File"}
        }
      }
    doc: "tumor BAM or CRAM"
  input_tumor_name: string
  input_normal_aligned:
    type: File
    secondaryFiles: |
      ${
        var dpath = self.location.replace(self.basename, "")
        if(self.nameext == '.bam'){
          return {"location": dpath+self.nameroot+".bai", "class": "File"}
        }
        else{
          return {"location": dpath+self.basename+".crai", "class": "File"}
        }
      }
    doc: "normal BAM or CRAM"
  input_normal_name: string
  manta_small_indels: { type: 'File?', doc: "Small indels output from Manta workflow", secondaryFiles: [.tbi] }
  use_manta_small_indels: {type: 'boolean?', default: false, doc: "Should the program use the small indels output from Manta in Strelka2 calling?"}
  vep_cache: { type: File, doc: "tar gzipped cache from ensembl/local converted cache", sbg:suggestedValue: {class: File, path: 607713829360f10e3982a425, name: homo_sapiens_vep_93_GRCh38.tar.gz} }
  hg38_strelka_bed: { type: File, doc: "Bgzipped interval bed file. Recommend padding 100bp for WXS; Recommend canonical chromosomes for WGS" , sbg:suggestedValue: {
      class: File, path: 5f500135e4b0370371c051ae, name: hg38_strelka.bed.gz} }
  hg38_strelka_tbi: { type: 'File?', doc: "Tabix index for hg38_strelka_bed", sbg:suggestedValue: {
      class: File, path: 5f500135e4b0370371c051aa, name: hg38_strelka.bed.gz.tbi} }
  output_basename: { type: string, doc: "String value to use as basename for outputs" }
  wgs_or_wxs: { type: { type: enum, name: sex, symbols: ["WGS", "WXS"] }, doc: "Select if this run is WGS or WXS" }

  # Optional with One Default
  select_vars_mode: { type: ['null', { type: enum, name: select_vars_mode, symbols: ["gatk", "grep"] }], default: "gatk", doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression" }
  vep_ref_build: { type: 'string?', default: "GRCh38", doc: "Genome ref build used, should line up with cache" }

  # Optional with Multiple Defaults (handled in choose_defaults)
  exome_flag: { type: 'string?', doc: "Whether to run in exome mode for callers. Y for WXS, N for WGS" }

  # annotation vars
  genomic_hotspots: { type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots", sbg:suggestedValue: [{class: File, path: 607713829360f10e3982a423, name: tert.bed}] }
  protein_snv_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots", sbg:suggestedValue: [{class: File, path: 607713829360f10e3982a426, name: protein_snv_cancer_hotspots_v2.tsv}] }
  protein_indel_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspots", sbg:suggestedValue: [{class: File, path: 607713829360f10e3982a424, name: protein_indel_cancer_hotspots_v2.tsv}] }
  retain_info: {type: string?, doc: "csv string with INFO fields that you want to keep", default: "MQ,MQ0,QSI,HotSpotAllele"}
  retain_fmt: {type: string?, doc: "csv string with FORMAT fields that you want to keep"}
  add_common_fields: {type: boolean?, doc: "Set to true if input is a strelka2 vcf that hasn't had common fields added", default: true}
  bcftools_annot_columns: {type: string, doc: "csv string of columns from annotation to port into the input vcf, i.e INFO/AF", default: "INFO/AF"}
  bcftools_annot_vcf: {type: File, doc: "bgzipped annotation vcf file", sbg:suggestedValue: {class: File, path: 5f50018fe4b054958bc8d2e3,
      name: af-only-gnomad.hg38.vcf.gz} }
  bcftools_annot_vcf_index: {type: File, doc: "index of bcftools_annot_vcf", sbg:suggestedValue: {class: File, path: 5f50018fe4b054958bc8d2e5,
      name: af-only-gnomad.hg38.vcf.gz.tbi}}
  bcftools_public_filter: {type: string?, doc: "Will hard filter final result to create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  gatk_filter_name: {type: 'string[]', doc: "Array of names for each filter tag to add, recommend: [\"NORM_DP_LOW\", \"GNOMAD_AF_HIGH\"]"}
  gatk_filter_expression: {type: 'string[]', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration, recommend: \"vc.getGenotype('\" + inputs.input_normal_name + \"').getDP() <= 7\"), \"AF > 0.001\"]"}
  disable_hotspot_annotation: { type: 'boolean?', doc: "Disable Hotspot Annotation and skip this task.", default: false }
  maf_center: {type: string?, doc: "Sequencing center of variant called", default: "."}

outputs:
  strelka2_prepass_vcf: { type: File, outputSource: run_strelka2/strelka2_prepass_vcf }
  strelka2_protected_outputs: { type: 'File[]', outputSource: run_strelka2/strelka2_protected_outputs }
  strelka2_public_outputs: { type: 'File[]', outputSource: run_strelka2/strelka2_public_outputs }

steps:
  choose_defaults:
    run: ../tools/mode_defaults.cwl
    in:
      input_mode: wgs_or_wxs
      exome_flag: exome_flag
    out: [out_exome_flag,out_cnvkit_wgs_mode,out_i_flag,out_lancet_padding,out_lancet_window,out_vardict_padding]

  prepare_reference:
    run: ../sub_workflows/prepare_reference.cwl
    in:
      input_fasta: reference_fasta
      input_fai: reference_fai
      input_dict: reference_dict
    out: [indexed_fasta,reference_dict]

  index_strelka_bed:
    run: ../tools/tabix_index.cwl
    in:
      input_file: hg38_strelka_bed
      input_index: hg38_strelka_tbi
    out: [output]

  index_bcftools_annot_vcf:
    run: ../tools/tabix_index.cwl
    in:
      input_file: bcftools_annot_vcf
      input_index: bcftools_annot_vcf_index
    out: [output]

  run_strelka2:
    run: ../sub_workflows/kfdrc_strelka2_sub_wf.cwl
    in:
      indexed_reference_fasta: prepare_reference/indexed_fasta
      reference_dict: prepare_reference/reference_dict
      hg38_strelka_bed: index_strelka_bed/output
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      manta_small_indels: manta_small_indels
      use_manta_small_indels: use_manta_small_indels
      exome_flag: exome_flag
      vep_cache: vep_cache
      vep_ref_build: vep_ref_build
      output_basename: output_basename
      select_vars_mode: select_vars_mode
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      retain_info: retain_info
      retain_fmt: retain_fmt
      add_common_fields: add_common_fields
      bcftools_annot_columns: bcftools_annot_columns
      bcftools_annot_vcf: index_bcftools_annot_vcf/output
      bcftools_public_filter: bcftools_public_filter
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      disable_hotspot_annotation: disable_hotspot_annotation
      maf_center: maf_center
    out:
      [strelka2_prepass_vcf, strelka2_protected_outputs, strelka2_public_outputs]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 2
