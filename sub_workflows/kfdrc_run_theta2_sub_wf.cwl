cwlVersion: v1.0
class: Workflow
id: kfdrc_cnvkit_theta2_adjust_wf
doc: "This workflow is normally run as a subworkflow and takes existing cnvkit results and adjusts them based on calculated tumor purity."

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  tumor_cns: {type: 'File', doc: "CNVkit output cns file"}
  reference_cnn: {type: 'File', doc: "CNVkit output cnn file"}
  tumor_sample_name: string
  normal_sample_name: string
  paired_vcf: {type: 'File', doc: "Combined somatic and germline call file. VarDict input recommended."}
  combined_include_expression: {type: ['null', string], doc: "Filter expression if vcf has non-PASS combined calls, use as-needed"}
  combined_exclude_expression: {type: ['null', string], doc: "Filter expression if vcf has non-PASS combined calls, use as-needed"}  
  min_theta2_frac: {type: ['null', float], doc: "Minimum fraction of genome with copy umber alterations.  Default is 0.05, recommend 0.01", default: 0.01}
  output_basename: string

outputs:
  theta2_adjusted_cns: {type: 'File?', outputSource: cnvkit_import_theta2/theta2_adjusted_cns}
  theta2_adjusted_seg: {type: 'File?', outputSource: cnvkit_import_theta2/theta2_adjusted_seg}
  theta2_subclonal_results: {type: ['null', 'File[]'], outputSource: [run_theta2/n3_graph, run_theta2/n2_results, run_theta2/best_results]}
  theta2_subclonal_cns: {type: ['null', 'File[]'], outputSource: cnvkit_import_theta2/theta2_subclone_cns}
  theta2_subclone_seg: {type: ['null', 'File[]'], outputSource: cnvkit_import_theta2/theta2_subclone_seg}

steps:
  bcftools_filter_combined_vcf:
    run: ../tools/bcftools_filter_vcf.cwl
    in:
      input_vcf: paired_vcf
      include_expression: combined_include_expression
      exclude_expression: combined_exclude_expression
      output_basename: output_basename
    out:
      [filtered_vcf]
  
  cnvkit_export_theta2:
    run: ../tools/cnvkit_export_theta2.cwl
    in:
      tumor_cns: tumor_cns
      reference_cnn: reference_cnn
      paired_vcf: bcftools_filter_combined_vcf/filtered_vcf
      normal_ID: normal_sample_name
      tumor_ID: tumor_sample_name
    out:
      [call_interval_count, call_tumor_snp, call_normal_snp]

  run_theta2:
    run: ../tools/theta2_purity.cwl
    in:
      tumor_snp: cnvkit_export_theta2/call_tumor_snp
      normal_snp: cnvkit_export_theta2/call_normal_snp
      interval_count: cnvkit_export_theta2/call_interval_count
      output_basename: output_basename
      min_frac: min_theta2_frac
    out:
      [n2_graph, n2_results, n2_withBounds, n3_graph, n3_results, n3_withBounds, best_results]
  
  cnvkit_import_theta2:
    run: ../tools/cnvkit_import_theta2.cwl
    in:
      tumor_cns: tumor_cns
      tumor_sample_name: tumor_sample_name
      output_basename: output_basename
      theta2_best_results: run_theta2/best_results
      theta2_n2_results: run_theta2/n2_results
    out:
      [theta2_adjusted_cns, theta2_adjusted_seg, theta2_subclone_cns, theta2_subclone_seg]
