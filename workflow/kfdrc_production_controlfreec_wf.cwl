cwlVersion: v1.0
class: Workflow
id: kfdrc_production_controlfreec_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
inputs:
  # Required
  reference_fasta: { type: File }
  reference_fai: { type: File }
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
  
  mate_copynumber_file_control: {type: File?, doc: "Normal cpn file from previous run. If used, will override bam use"}
  mate_copynumber_file_sample: {type: File?, doc: "Tumor cpn file from previous run. If used, will override bam use"}
  gem_mappability_file: {type: File?, doc: "GEM mappability file to make read count adjustments with"}
  min_subclone_presence: {type: float?, doc: "Use if you want to detect sublones. Recommend 0.2 for WGS, 0.3 for WXS"}
  cfree_chr_len: { type: File, doc: "file with chromosome lengths" }
  cfree_ploidy: { type: 'int[]', doc: "Array of ploidy possibilities for ControlFreeC to try" }
  output_basename: { type: string, doc: "String value to use as basename for outputs" }
  wgs_or_wxs: { type: { type: enum, name: sex, symbols: ["WGS", "WXS"] }, doc: "Select if this run is WGS or WXS" }

  # Optional with One Default
  cfree_threads: { type: 'int?', default: 16, doc: "For ControlFreeC. Recommend 16 max, as I/O gets saturated after that losing any advantage" }
  cfree_mate_orientation_control: { type: ['null', { type: enum, name: mate_orientation_control, symbols: ["0", "FR", "RF", "FF"] }], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)" }
  cfree_mate_orientation_sample: { type: ['null', { type: enum, name: mate_orientation_sample, symbols: ["0", "FR", "RF", "FF"] }], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)" }

  # Optional with Multiple Defaults (handled in choose_defaults)
  i_flag: { type: 'string?', doc: "Flag to intersect germline calls on padded regions. Use N if you want to skip this or have a WGS run" }

  # Optional
  b_allele: { type: 'File?', doc: "germline calls, needed for BAF.  GATK HC VQSR input recommended.  Tool will prefilter for germline and pass if expression given" }
  b_allele_index: { type: 'File?', doc: "Tabix index for b_allele" }
  cfree_coeff_var: { type: 'float?', default: 0.05, doc: "Coefficient of variation to set window size.  Default 0.05 recommended" }
  cfree_contamination_adjustment: { type: 'boolean?', doc: "TRUE or FALSE to have ControlFreec estimate normal contam" }
  cfree_sex: { type: ['null', { type: enum, name: sex, symbols: ["XX", "XY"] }], doc: "If known, XX for female, XY for male" }

  # WXS only Fields
  unpadded_capture_regions: { type: 'File?', doc: "Capture regions with NO padding for cnv calling" }

outputs:
  ctrlfreec_cnvs: {type: File, outputSource: run_controlfreec/ctrlfreec_cnvs}
  ctrlfreec_pval: { type: File, outputSource: run_controlfreec/ctrlfreec_pval }
  ctrlfreec_config: { type: File, outputSource: run_controlfreec/ctrlfreec_config }
  ctrlfreec_pngs: { type: 'File[]', outputSource: run_controlfreec/ctrlfreec_pngs }
  ctrlfreec_bam_ratio: { type: File, outputSource: run_controlfreec/ctrlfreec_bam_ratio }
  ctrlfreec_bam_seg: { type: File, outputSource: run_controlfreec/ctrlfreec_bam_seg }
  ctrlfreec_baf: { type: File, outputSource: run_controlfreec/ctrlfreec_baf }
  ctrlfreec_info: { type: File, outputSource: run_controlfreec/ctrlfreec_info }

steps:
  choose_defaults:
    run: ../tools/mode_defaults.cwl
    in:
      input_mode: wgs_or_wxs
      i_flag: i_flag
    out: [out_exome_flag,out_cnvkit_wgs_mode,out_i_flag,out_lancet_padding,out_lancet_window,out_vardict_padding]

  prepare_reference:
    run: ../sub_workflows/prepare_reference.cwl
    in:
      input_fasta: reference_fasta
      input_fai: reference_fai
      input_dict: reference_dict
    out: [indexed_fasta,reference_dict]

  index_b_allele:
    run: ../tools/tabix_index.cwl
    in:
      input_file: b_allele
      input_index: b_allele_index
    out: [output]

  bedtools_intersect_germline:
    run: ../tools/bedtools_intersect.cwl
    in:
      input_vcf: index_b_allele/output
      output_basename: output_basename
      input_bed_file: unpadded_capture_regions
      flag: choose_defaults/out_i_flag
    out:
      [intersected_vcf]

  gatk_filter_germline:
    run: ../tools/gatk_filter_germline_variant.cwl
    in:
      input_vcf: bedtools_intersect_germline/intersected_vcf
      reference_fasta: prepare_reference/indexed_fasta
      output_basename: output_basename
    out:
      [filtered_vcf, filtered_pass_vcf]

  samtools_cram2bam_plus_calmd_tumor:
    run: ../tools/samtools_cram2bam_plus_calmd.cwl
    in:
      input_reads: input_tumor_aligned
      threads:
        valueFrom: ${return 16;}
      reference: prepare_reference/indexed_fasta
    out: [bam_file]

  samtools_cram2bam_plus_calmd_normal:
    run: ../tools/samtools_cram2bam_plus_calmd.cwl
    in:
      input_reads: input_normal_aligned
      threads:
        valueFrom: ${return 16;}
      reference: prepare_reference/indexed_fasta
    out: [bam_file]

  run_controlfreec:
    run: ../sub_workflows/kfdrc_controlfreec_sub_wf.cwl
    in:
      mate_copynumber_file_control: mate_copynumber_file_control
      mate_copynumber_file_sample: mate_copynumber_file_sample
      gem_mappability_file: gem_mappability_file
      min_subclone_presence: min_subclone_presence
      input_tumor_aligned: samtools_cram2bam_plus_calmd_tumor/bam_file
      input_tumor_name: input_tumor_name
      input_normal_aligned: samtools_cram2bam_plus_calmd_normal/bam_file
      threads: cfree_threads
      output_basename: output_basename
      ploidy: cfree_ploidy
      mate_orientation_sample: cfree_mate_orientation_sample
      mate_orientation_control: cfree_mate_orientation_control
      capture_regions: unpadded_capture_regions
      indexed_reference_fasta: prepare_reference/indexed_fasta
      reference_fai: reference_fai
      b_allele: gatk_filter_germline/filtered_pass_vcf
      chr_len: cfree_chr_len
      coeff_var: cfree_coeff_var
      contamination_adjustment: cfree_contamination_adjustment
      cfree_sex: cfree_sex
    out:
      [ctrlfreec_cnvs, ctrlfreec_pval, ctrlfreec_config, ctrlfreec_pngs, ctrlfreec_bam_ratio, ctrlfreec_bam_seg, ctrlfreec_baf, ctrlfreec_info]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4
