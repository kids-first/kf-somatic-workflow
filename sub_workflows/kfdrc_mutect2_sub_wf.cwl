cwlVersion: v1.2
class: Workflow
id: kfdrc_mutect2_sub_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  indexed_reference_fasta: {type: 'File', secondaryFiles: [.fai, ^.dict]}
  reference_dict: File
  bed_invtl_split: {type: 'File[]', doc: "Bed file intervals passed on from and outside pre-processing step"}
  af_only_gnomad_vcf: {type: 'File', secondaryFiles: ['.tbi']}
  exac_common_vcf: {type: 'File', secondaryFiles: ['.tbi']}
  input_tumor_aligned: { type: File, secondaryFiles: [{pattern: ".bai", required: false},
      {pattern: "^.bai", required: false}, {pattern: ".crai", required: false}, {
        pattern: "^.crai", required: false}], doc: "tumor BAM or CRAM" }
  input_tumor_name: string
  old_tumor_name: { type: 'string?', doc: "If `SM:` sample name in te align file is different than `input_tumor_name`, you **must** provide it here"}
  input_normal_aligned: { type: File, secondaryFiles: [{pattern: ".bai", required: false},
      {pattern: "^.bai", required: false}, {pattern: ".crai", required: false}, {
        pattern: "^.crai", required: false}], doc: "normal BAM or CRAM" }
  input_normal_name: string
  old_normal_name: { type: 'string?', doc: "If `SM:` sample name in te align file is different than `input_normal_name`, you **must** provide it here"}
  exome_flag: {type: ['null', string], doc: "set to 'Y' for exome mode"}
  select_vars_mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression", default: "gatk"}
  tool_name: {type: 'string?', doc: "String to describe what tool was run as part of file name", default: "mutect2_somatic"}
  output_basename: string
  getpileup_memory: {type: 'int?'}
  learnorientation_memory: {type: 'int?'}
  filtermutectcalls_memory: {type: 'int?'}
  vep_cache: {type: 'File', doc: "tar gzipped cache from ensembl/local converted cache"}
  vep_ram: {type: 'int?', doc: "In GB, may need to increase this value depending on the size/complexity of input"}
  vep_cores: {type: 'int?', doc: "Number of cores to use. May need to increase for really large inputs"}
  vep_buffer_size: {type: 'int?', doc: "Increase or decrease to balance speed and memory usage"}
  dbnsfp: { type: 'File?', secondaryFiles: [.tbi,^.readme.txt], doc: "VEP-formatted plugin file, index, and readme file containing dbNSFP annotations" }
  dbnsfp_fields: { type: 'string?', doc: "csv string with desired fields to annotate. Use ALL to grab all"}
  merged: { type: 'boolean?', doc: "Set to true if merged cache used", default: true }
  cadd_indels: { type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD indel annotations" }
  cadd_snvs: { type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD SNV annotations" }
  run_cache_existing: { type: boolean, doc: "Run the check_existing flag for cache" }
  run_cache_af: { type: boolean, doc: "Run the allele frequency flags for cache" }

  # annotation vars
  genomic_hotspots: { type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots" }
  protein_snv_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots" }
  protein_indel_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspots" }
  retain_info: {type: 'string?', doc: "csv string with INFO fields that you want to keep", default: "MBQ,TLOD,HotSpotAllele"}
  retain_fmt: {type: 'string?', doc: "csv string with FORMAT fields that you want to keep"}
  retain_ann: { type: 'string?', doc: "csv string of annotations (within the VEP CSQ/ANN) to retain as extra columns in MAF" }
  add_common_fields: {type: 'boolean?', doc: "Set to true if input is a strelka2 vcf that hasn't had common fields added", default: false}
  bcftools_annot_columns: {type: 'string?', doc: "csv string of columns from annotation to port into the input vcf, i.e INFO/AF", default: "INFO/AF"}
  bcftools_strip_columns: {type: 'string?', doc: "csv string of columns to strip if needed to avoid conflict, i.e INFO/AF"}
  bcftools_annot_vcf: {type: 'File?', secondaryFiles: ['.tbi'], doc: "additional bgzipped annotation vcf file"}
  bcftools_public_filter: {type: 'string?', doc: "Will hard filter final result to create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  gatk_filter_name: {type: 'string[]', doc: "Array of names for each filter tag to add, recommend: [\"NORM_DP_LOW\", \"GNOMAD_AF_HIGH\"]"}
  gatk_filter_expression: {type: 'string[]', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration, recommend: \"vc.getGenotype('\" + inputs.input_normal_name + \"').getDP() <= 7\"), \"AF > 0.001\"]"}
  disable_hotspot_annotation: { type: 'boolean?', doc: "Disable Hotspot Annotation and skip this task.", default: false }
  maf_center: {type: 'string?', doc: "Sequencing center of variant called", default: "."}

outputs:
  mutect2_filtered_stats: {type: 'File', outputSource: filter_mutect2_vcf/stats_table}
  mutect2_filtered_vcf: { type: 'File', outputSource: pickvalue_workaround/output }
  mutect2_protected_outputs: {type: 'File[]', outputSource: annotate/annotated_protected}
  mutect2_public_outputs: {type: 'File[]', outputSource: annotate/annotated_public}

steps:
  mutect2:
    run: ../tools/gatk_Mutect2.cwl
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge
    in:
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: 
        source: [old_tumor_name, input_tumor_name]
        pickValue: first_non_null
      input_normal_aligned: input_normal_aligned
      input_normal_name:
        source: [old_normal_name, input_normal_name]
        pickValue: first_non_null
      reference: indexed_reference_fasta
      interval_list: bed_invtl_split
      af_only_gnomad_vcf: af_only_gnomad_vcf
      exome_flag: exome_flag
    scatter: [interval_list]
    out: [mutect2_vcf, f1r2_counts, mutect_stats]

  mutect2_filter_support:
    run: ../sub_workflows/kfdrc_mutect2_filter_support_subwf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      reference_dict: reference_dict
      wgs_calling_interval_list: bed_invtl_split
      input_tumor_aligned: input_tumor_aligned
      input_normal_aligned: input_normal_aligned
      exac_common_vcf: exac_common_vcf
      output_basename: output_basename
      f1r2_counts: mutect2/f1r2_counts
      getpileup_memory: getpileup_memory
      learnorientation_memory: learnorientation_memory
    out: [contamination_table, segmentation_table, f1r2_bias]

  merge_mutect2_vcf:
    run: ../tools/gatk_mergevcfs.cwl
    label: Merge mutect2 vcf
    in:
      input_vcfs: mutect2/mutect2_vcf
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name: tool_name
    out: [merged_vcf]

  merge_mutect2_stats:
    run: ../tools/gatk_mergemutectstats.cwl
    label: Merge mutect2 stats
    in:
      input_stats: mutect2/mutect_stats
      output_basename: output_basename
    out: [merged_stats]

  filter_mutect2_vcf:
    run: ../tools/gatk_filtermutectcalls.cwl
    in:
      mutect_vcf: merge_mutect2_vcf/merged_vcf
      mutect_stats: merge_mutect2_stats/merged_stats
      reference: indexed_reference_fasta
      output_basename: output_basename
      contamination_table: mutect2_filter_support/contamination_table
      segmentation_table: mutect2_filter_support/segmentation_table
      ob_priors: mutect2_filter_support/f1r2_bias
      max_memory: filtermutectcalls_memory
    out: [stats_table, filtered_vcf]

  rename_vcf_samples:
    when: $(inputs.old_tumor_name != null)
    run: ../tools/bcftools_reheader_vcf.cwl
    in:
      input_vcf: filter_mutect2_vcf/filtered_vcf
      input_normal_name: input_normal_name
      input_tumor_name: input_tumor_name
      old_tumor_name: old_tumor_name
    out: [reheadered_vcf]

  pickvalue_workaround:
    run: ../tools/expression_pickvalue_workaround.cwl
    in:
      input_file:
        source: [rename_vcf_samples/reheadered_vcf, filter_mutect2_vcf/filtered_vcf]
        pickValue: first_non_null
    out: [output]

  gatk_selectvariants_mutect2:
    run: ../tools/gatk_selectvariants.cwl
    label: GATK Select PASS
    in:
      input_vcf: pickvalue_workaround/output
      output_basename: output_basename
      tool_name: tool_name
      mode: select_vars_mode
    out: [pass_vcf]

  annotate:
    run: ../workflow/kfdrc_annot_vcf_wf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf: gatk_selectvariants_mutect2/pass_vcf
      input_tumor_name: input_tumor_name
      input_normal_name: input_normal_name
      add_common_fields: add_common_fields
      retain_info: retain_info
      retain_fmt: retain_fmt
      retain_ann: retain_ann
      bcftools_annot_columns: bcftools_annot_columns
      bcftools_strip_columns: bcftools_strip_columns
      bcftools_annot_vcf: bcftools_annot_vcf
      bcftools_public_filter: bcftools_public_filter
      dbnsfp: dbnsfp
      dbnsfp_fields: dbnsfp_fields
      merged: merged
      cadd_indels: cadd_indels
      cadd_snvs: cadd_snvs
      run_cache_af: run_cache_af
      run_cache_existing: run_cache_existing
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      vep_cache: vep_cache
      vep_ram: vep_ram
      vep_cores: vep_cores
      vep_buffer_size: vep_buffer_size
      disable_hotspot_annotation: disable_hotspot_annotation
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      maf_center: maf_center
      output_basename: output_basename
      tool_name: tool_name
    out: [annotated_protected, annotated_public]

$namespaces:
  sbg: https://sevenbridges.com
