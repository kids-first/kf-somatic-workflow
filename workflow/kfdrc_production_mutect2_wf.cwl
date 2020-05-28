cwlVersion: v1.0
class: Workflow
id: kfdrc_production_mutect2_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
inputs:
  # Required
  reference_fasta: {type: File }
  reference_fai: { type: 'File?' }
  reference_dict: { type: 'File?' }
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
  vep_cache: {type: File, doc: "tar gzipped cache from ensembl/local converted cache"}
  mutect2_af_only_gnomad_vcf: {type: File}
  mutect2_af_only_gnomad_tbi: {type: 'File?', doc: "Tabix index for mutect2_af_only_gnomad_vcf"}
  mutect2_exac_common_vcf: {type: File}
  mutect2_exac_common_tbi: {type: 'File?', doc: "Tabix index for mutect2_exac_common_vcf"}
  output_basename: {type: string, doc: "String value to use as basename for outputs"}
  wgs_or_wxs: {type: {type: enum, name: sex, symbols: ["WGS", "WXS"] }, doc: "Select if this run is WGS or WXS"}

  # Optional with One Default
  select_vars_mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], default: "gatk", doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression"}
  vep_ref_build: {type: 'string?', default: "GRCh38", doc: "Genome ref build used, should line up with cache"}

  # Optional with Multiple Defaults (handled in choose_defaults)
  exome_flag: {type: string?, doc: "Whether to run in exome mode for callers. Y for WXS, N for WGS"}

  # WGS only Fields
  wgs_calling_interval_list: {type: File?, doc: "GATK intervals list-style, or bed file.  Recommend canocical chromosomes with N regions removed"}

  # WXS only Fields
  padded_capture_regions: {type: 'File?', doc: "Recommend 100bp pad, for somatic variant"}

outputs:
  mutect2_vep_vcf: {type: File, outputSource: run_mutect2/mutect2_vep_vcf}
  mutect2_vep_tbi: {type: File, outputSource: run_mutect2/mutect2_vep_tbi}
  mutect2_prepass_vcf: {type: File, outputSource: run_mutect2/mutect2_filtered_vcf}
  mutect2_vep_maf: {type: File, outputSource: run_mutect2/mutect2_vep_maf}

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

  index_mutect_gnomad:
    run: ../tools/tabix_index.cwl
    in:
      input_file: mutect2_af_only_gnomad_vcf
      input_index: mutect2_af_only_gnomad_tbi
    out: [output]

  index_mutect_exac:
    run: ../tools/tabix_index.cwl
    in:
      input_file: mutect2_exac_common_vcf
      input_index: mutect2_exac_common_tbi
    out: [output]

  select_interval_list:
    run: ../tools/mode_selector.cwl
    in:
      input_mode: wgs_or_wxs
      wgs_input: wgs_calling_interval_list
      wxs_input: padded_capture_regions
    out: [output]

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

  select_mutect_bed_interval:
    run: ../tools/mode_selector.cwl
    in:
      input_mode: wgs_or_wxs
      wgs_input: gatk_intervallisttools/output
      wxs_input: gatk_intervallisttools/output
    out: [output]

  run_mutect2:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge
    run: ../sub_workflows/kfdrc_mutect2_sub_wf.cwl
    in:
      indexed_reference_fasta: prepare_reference/indexed_fasta
      reference_dict: prepare_reference/reference_dict
      bed_invtl_split: select_mutect_bed_interval/output
      af_only_gnomad_vcf: index_mutect_gnomad/output
      exac_common_vcf: index_mutect_exac/output
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      exome_flag: choose_defaults/out_exome_flag
      vep_cache: vep_cache
      vep_ref_build: vep_ref_build
      output_basename: output_basename
      select_vars_mode: select_vars_mode
    out:
      [mutect2_filtered_stats, mutect2_filtered_vcf, mutect2_vep_vcf, mutect2_vep_tbi, mutect2_vep_maf]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 6
  - class: 'sbg:AWSInstanceType'
    value: c5.9xlarge
