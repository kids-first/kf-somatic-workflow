cwlVersion: v1.0
class: Workflow
id: kfdrc_production_vardict_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
inputs:
  # Required
  reference_fasta: {type: 'File', "sbg:suggestedValue": {class: File, path: 60639014357c3a53540ca7a3,
      name: Homo_sapiens_assembly38.fasta}}
  reference_fai: {type: 'File?', "sbg:suggestedValue": {class: File, path: 60639016357c3a53540ca7af,
      name: Homo_sapiens_assembly38.fasta.fai}}
  reference_dict: {type: 'File?', "sbg:suggestedValue": {class: File, path: 60639019357c3a53540ca7e7,
      name: Homo_sapiens_assembly38.dict}}
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
  output_basename: { type: 'string', doc: "String value to use as basename for outputs" }
  wgs_or_wxs: { type: { type: enum, name: sex, symbols: ["WGS", "WXS"] }, doc: "Select if this run is WGS or WXS" }

  # Optional with One Default
  select_vars_mode: { type: ['null', { type: enum, name: select_vars_mode, symbols: ["gatk", "grep"] }], default: "gatk", doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression" }
  vardict_cpus: { type: 'int?', default: 9, doc: "Number of CPUs for Vardict to use" }
  vardict_min_vaf: { type: 'float?', default: 0.05, doc: "Min variant allele frequency for vardict to consider. Recommend 0.05" }
  vardict_ram: { type: 'int?', default: 18, doc: "GB of RAM to allocate to Vardict" }

  # Optional with Multiple Defaults (handled in choose_defaults)
  exome_flag: { type: 'string?', doc: "Whether to run in exome mode for callers. Y for WXS, N for WGS" }
  vardict_padding: { type: 'int?', doc: "Padding to add to input intervals, recommend 0 if intervals already padded such as in WXS, 150 if not such as in WGS" }

  # WGS only Fields
  wgs_calling_interval_list: {type: 'File?', doc: "GATK intervals list-style, or bed\
      \ file.  Recommend canocical chromosomes with N regions removed", "sbg:suggestedValue": {
      class: File, path: 5f500135e4b0370371c051b6, name: wgs_canonical_calling_regions.hg38.bed}}

  # WXS only Fields
  padded_capture_regions: { type: 'File?', doc: "Recommend 100bp pad, for somatic variant" }

  # VEP params
  vep_cache: {type: 'File', doc: "tar gzipped cache from ensembl/local converted cache",  "sbg:suggestedValue": {class: File, path: 6332f8e47535110eb79c794f,
      name: homo_sapiens_merged_vep_105_indexed_GRCh38.tar.gz}}
  vep_ram: {type: 'int?', default: 32, doc: "In GB, may need to increase this value depending on the size/complexity of input"}
  vep_cores: {type: 'int?', default: 16, doc: "Number of cores to use. May need to increase for really large inputs"}
  vep_buffer_size: {type: 'int?', default: 1000, doc: "Increase or decrease to balance speed and memory usage"}
  dbnsfp: { type: 'File?', secondaryFiles: [.tbi,^.readme.txt], doc: "VEP-formatted plugin file, index, and readme file containing dbNSFP annotations" }
  dbnsfp_fields: { type: 'string?', doc: "csv string with desired fields to annotate. Use ALL to grab all"}
  merged: { type: 'boolean?', doc: "Set to true if merged cache used", default: true }
  cadd_indels: { type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD indel annotations" }
  cadd_snvs: { type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD SNV annotations" }
  run_cache_existing: { type: 'boolean?', doc: "Run the check_existing flag for cache" }
  run_cache_af: { type: 'boolean?', doc: "Run the allele frequency flags for cache" }

  # annotation vars
  genomic_hotspots: { type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots", "sbg:suggestedValue": [{class: File, path: 607713829360f10e3982a423, name: tert.bed}] }
  protein_snv_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots", "sbg:suggestedValue": [{class: File, path: 607713829360f10e3982a426, name: protein_snv_cancer_hotspots_v2.tsv}] }
  protein_indel_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspots", "sbg:suggestedValue": [{class: File, path: 607713829360f10e3982a424, name: protein_indel_cancer_hotspots_v2.tsv}] }
  retain_info: {type: 'string?', doc: "csv string with INFO fields that you want to keep", default: "gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,MSI,MSILEN,SOR,SSF,HotSpotAllele"}
  retain_fmt: {type: 'string?', doc: "csv string with FORMAT fields that you want to keep"}
  retain_ann: { type: 'string?', doc: "csv string of annotations (within the VEP CSQ/ANN) to retain as extra columns in MAF", default: "HGVSg" }
  bcftools_annot_columns: {type: 'string?', doc: "csv string of columns from annotation to port into the input vcf, i.e INFO/AF", default: "INFO/AF"}
  bcftools_strip_columns: {type: 'string?', doc: "csv string of columns to strip if needed to avoid conflict, i.e INFO/AF"}
  bcftools_annot_vcf: {type: 'File', doc: "bgzipped annotation vcf file", "sbg:suggestedValue": {
      class: File, path: 6324ef5ad01163633daa00d8, name: gnomad_3.1.1.vwb_subset.vcf.gz}}
  bcftools_annot_vcf_index: {type: 'File', doc: "index of bcftools_annot_vcf", "sbg:suggestedValue": {
      class: File, path: 6324ef5ad01163633daa00d7, name: gnomad_3.1.1.vwb_subset.vcf.gz.tbi}}
  bcftools_public_filter: {type: 'string?', doc: "Will hard filter final result to create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  gatk_filter_name: {type: 'string[]', doc: "Array of names for each filter tag to add, recommend: [\"NORM_DP_LOW\", \"GNOMAD_AF_HIGH\"]"}
  gatk_filter_expression: {type: 'string[]', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration, recommend: \"vc.getGenotype('\" + inputs.input_normal_name + \"').getDP() <= 7\"), \"AF > 0.001\"]"}
  disable_hotspot_annotation: { type: 'boolean?', doc: "Disable Hotspot Annotation and skip this task.", default: false }
  maf_center: {type: 'string?', doc: "Sequencing center of variant called", default: "."}

outputs:
  vardict_prepass_vcf: { type: 'File', outputSource: run_vardict/vardict_prepass_vcf }
  vardict_protected_outputs: { type: 'File[]', outputSource: run_vardict/vardict_protected_outputs }
  vardict_public_outputs: { type: 'File[]', outputSource: run_vardict/vardict_public_outputs }

steps:
  choose_defaults:
    run: ../tools/mode_defaults.cwl
    in:
      input_mode: wgs_or_wxs
      exome_flag: exome_flag
      vardict_padding: vardict_padding
    out: [out_exome_flag,out_cnvkit_wgs_mode,out_i_flag,out_lancet_padding,out_lancet_window,out_vardict_padding]

  prepare_reference:
    run: ../sub_workflows/prepare_reference.cwl
    in:
      input_fasta: reference_fasta
      input_fai: reference_fai
      input_dict: reference_dict
    out: [indexed_fasta,reference_dict]

  index_bcftools_annot_vcf:
    run: ../tools/tabix_index.cwl
    in:
      input_file: bcftools_annot_vcf
      input_index: bcftools_annot_vcf_index
    out: [output]

  select_interval_list:
    run: ../tools/mode_selector.cwl
    in:
      input_mode: wgs_or_wxs
      wgs_input: wgs_calling_interval_list
      wxs_input: padded_capture_regions
    out: [output]

  python_vardict_interval_split:
    run: ../tools/python_vardict_interval_split.cwl
    doc: "Custom interval list generation for vardict input. Briefly, ~60M bp per interval list, 20K bp intervals, lists break on chr and N reginos only"
    in:
      wgs_bed_file: select_interval_list/output
    out: [split_intervals_bed]

  gatk_intervallisttools:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: select_interval_list/output
      reference_dict: prepare_reference/reference_dict
      exome_flag: choose_defaults/out_exome_flag
      scatter_ct:
        valueFrom: ${return 50}
      bands:
        valueFrom: ${return 80000000}
    out: [output]

  samtools_cram2bam_plus_calmd_tumor:
    run: ../tools/samtools_calmd.cwl
    in:
      input_reads: input_tumor_aligned
      reference: prepare_reference/indexed_fasta
    out: [bam_file]

  samtools_cram2bam_plus_calmd_normal:
    run: ../tools/samtools_calmd.cwl
    in:
      input_reads: input_normal_aligned
      reference: prepare_reference/indexed_fasta
    out: [bam_file]

  select_vardict_bed_interval:
    run: ../tools/mode_list_selector.cwl
    in:
      input_mode: wgs_or_wxs
      wgs_input: python_vardict_interval_split/split_intervals_bed
      wxs_input: gatk_intervallisttools/output
    out: [output]

  run_vardict:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge
    run: ../sub_workflows/kfdrc_vardict_sub_wf.cwl
    in:
      indexed_reference_fasta: prepare_reference/indexed_fasta
      input_tumor_aligned: samtools_cram2bam_plus_calmd_tumor/bam_file
      input_tumor_name: input_tumor_name
      input_normal_aligned: samtools_cram2bam_plus_calmd_normal/bam_file
      input_normal_name: input_normal_name
      output_basename: output_basename
      reference_dict: prepare_reference/reference_dict
      bed_invtl_split: select_vardict_bed_interval/output
      padding: choose_defaults/out_vardict_padding
      min_vaf: vardict_min_vaf
      select_vars_mode: select_vars_mode
      cpus: vardict_cpus
      ram: vardict_ram
      retain_info: retain_info
      retain_fmt: retain_fmt
      retain_ann: retain_ann
      bcftools_annot_columns: bcftools_annot_columns
      bcftools_strip_columns: bcftools_strip_columns
      bcftools_annot_vcf: index_bcftools_annot_vcf/output
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
    out:
      [vardict_prepass_vcf, vardict_protected_outputs, vardict_public_outputs]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4
