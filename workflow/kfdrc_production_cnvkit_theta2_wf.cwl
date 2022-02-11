cwlVersion: v1.2
class: Workflow
id: kfdrc_production_cnvkit_theta2_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
inputs:
  # Required
  reference_fasta: { type: File }
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
  cnvkit_annotation_file: { type: File, doc: "refFlat.txt file" }
  cnvkit_cnn: {type: ['null', File], doc: "If running using an existing .cnn, supply here"}
  output_basename: { type: string, doc: "String value to use as basename for outputs" }
  wgs_or_wxs: { type: { type: enum, name: sex, symbols: ["WGS", "WXS"] }, doc: "Select if this run is WGS or WXS" }

  # Optional with Multiple Defaults (handled in choose_defaults)
  cnvkit_wgs_mode: { type: 'string?', doc: "for WGS mode, input Y. leave blank for WXS/hybrid mode" }
  i_flag: { type: 'string?', doc: "Flag to intersect germline calls on padded regions. Use N if you want to skip this or have a WGS run" }

  # Optional
  b_allele: { type: 'File?', doc: "germline calls, needed for BAF.  GATK HC VQSR input recommended.  Tool will prefilter for germline and pass if expression given" }
  b_allele_index: { type: 'File?', doc: "Tabix index for b_allele" }
  cnvkit_sex: {type: ['null', {type: enum, name: cnvkit_sex, symbols: ["x", "y"]}], doc: "Sex, for simplicity x for female y for male", default: "x"}

  # WXS only Fields
  unpadded_capture_regions: { type: 'File?', doc: "Capture regions with NO padding for cnv calling" }

  # theta inputs
  min_theta2_frac: { type: 'float?', default: 0.01, doc: "Minimum fraction of genome with copy umber alterations.  Default is 0.05, recommend 0.01" }
  combined_include_expression: { type: 'string?', doc: "Theta2 Purity value: Filter expression if vcf has non-PASS combined calls, use as-needed, i.e. for VarDict: FILTER=\"PASS\" && (INFO/STATUS=\"Germline\" | INFO/STATUS=\"StrongSomatic\")" }
  combined_exclude_expression: { type: 'string?', doc: "Theta2 Purity value: Filter expression if vcf has non-PASS combined calls, use as-needed" }
  vardict_prepass_vcf: { type: File, doc: "Prepass VCF output from Vardict" }


outputs:
  cnvkit_cnr: { type: File, outputSource: run_cnvkit/cnvkit_cnr }
  cnvkit_cnn_output: { type: ['null', File], outputSource: run_cnvkit/cnvkit_cnn_output }
  cnvkit_calls: { type: File, outputSource: run_cnvkit/cnvkit_calls }
  cnvkit_metrics: { type: File, outputSource: run_cnvkit/cnvkit_metrics }
  cnvkit_gainloss: { type: File, outputSource: run_cnvkit/cnvkit_gainloss }
  cnvkit_seg: { type: File, outputSource: run_cnvkit/cnvkit_seg }
  cnvkit_scatter_plot: {type: File, outputSource: run_cnvkit/cnvkit_scatter_plot}
  cnvkit_diagram: {type: File, outputSource: run_cnvkit/cnvkit_diagram}
  theta2_calls: { type: 'File?', outputSource: run_theta2_purity/theta2_adjusted_cns }
  theta2_seg: { type: 'File?', outputSource: run_theta2_purity/theta2_adjusted_seg }
  theta2_subclonal_results: {type: ['null', 'File[]'], outputSource: expression_flatten_subclonal_results/output}
  theta2_subclonal_cns: { type: ['null', 'File[]'], outputSource: run_theta2_purity/theta2_subclonal_cns }
  theta2_subclone_seg: { type: ['null', 'File[]'], outputSource: run_theta2_purity/theta2_subclone_seg }


steps:
  choose_defaults:
    run: ../tools/mode_defaults.cwl
    in:
      input_mode: wgs_or_wxs
      cnvkit_wgs_mode: cnvkit_wgs_mode
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
    run: ../tools/samtools_calmd.cwl
    in:
      input_reads: input_tumor_aligned
      threads:
        valueFrom: ${return 16;}
      reference: prepare_reference/indexed_fasta
    out: [bam_file]

  samtools_cram2bam_plus_calmd_normal:
    run: ../tools/samtools_calmd.cwl
    in:
      input_reads: input_normal_aligned
      threads:
        valueFrom: ${return 16;}
      reference: prepare_reference/indexed_fasta
    out: [bam_file]

  run_cnvkit:
    run: ../sub_workflows/kfdrc_cnvkit_sub_wf.cwl
    in:
      input_tumor_aligned: samtools_cram2bam_plus_calmd_tumor/bam_file
      tumor_sample_name: input_tumor_name
      input_normal_aligned: samtools_cram2bam_plus_calmd_normal/bam_file
      cnvkit_cnn_input: cnvkit_cnn
      reference: prepare_reference/indexed_fasta
      normal_sample_name: input_normal_name
      capture_regions: unpadded_capture_regions
      wgs_mode: choose_defaults/out_cnvkit_wgs_mode
      b_allele_vcf: gatk_filter_germline/filtered_pass_vcf
      annotation_file: cnvkit_annotation_file
      output_basename: output_basename
      sex: cnvkit_sex
    out:
      [cnvkit_cnr, cnvkit_cnn_output, cnvkit_calls, cnvkit_metrics, cnvkit_gainloss, cnvkit_seg, cnvkit_scatter_plot, cnvkit_diagram]

  run_theta2_purity:
    run: ../sub_workflows/kfdrc_run_theta2_sub_wf.cwl
    in:
      tumor_cns: run_cnvkit/cnvkit_calls
      reference_cnn:
        source: [cnvkit_cnn, run_cnvkit/cnvkit_cnn_output]
        pickValue: first_non_null
      tumor_sample_name: input_tumor_name
      normal_sample_name: input_normal_name
      paired_vcf: vardict_prepass_vcf
      combined_include_expression: combined_include_expression
      combined_exclude_expression: combined_exclude_expression
      min_theta2_frac: min_theta2_frac
      output_basename: output_basename
    out:
      [theta2_adjusted_cns, theta2_adjusted_seg, theta2_subclonal_results, theta2_subclonal_cns, theta2_subclone_seg]

  expression_flatten_subclonal_results:
    run: ../tools/expression_flatten_list.cwl
    in:
      input_list: run_theta2_purity/theta2_subclonal_results
    out: [output]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4
