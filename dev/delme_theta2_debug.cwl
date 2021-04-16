cwlVersion: v1.0
class: Workflow
id: delme-debug-workflow

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  tumor_cns: {type: File, doc: "CNVkit output cns file"}
  reference_cnn: {type: File, doc: "CNVkit output cnn file"}
  tumor_sample_name: string
  normal_sample_name: string
  paired_vcf: {type: File, doc: "Combined somatic and germline call file. VarDict input recommended."}
  combined_include_expression: {type: ['null', string], doc: "Filter expression if vcf has non-PASS combined calls, use as-needed"}
  combined_exclude_expression: {type: ['null', string], doc: "Filter expression if vcf has non-PASS combined calls, use as-needed"}  
  min_theta2_frac: {type: ['null', float], doc: "Minimum fraction of genome with copy umber alterations.  Default is 0.05, recommend 0.01", default: 0.01}
  output_basename: string

outputs:
  theta2_subclonal_results: {type: 'File[]?', outputSource: expression_flatten_subclobal_results/output}

steps:
  run_theta2_purity:
    run: ../sub_workflows/kfdrc_run_theta2_sub_wf.cwl
    in:
      tumor_cns: tumor_cns
      reference_cnn: reference_cnn
      tumor_sample_name: tumor_sample_name
      normal_sample_name: normal_sample_name
      paired_vcf: paired_vcf
      combined_include_expression: combined_include_expression
      combined_exclude_expression: combined_exclude_expression
      min_theta2_frac: min_theta2_frac
      output_basename: output_basename
    out: [theta2_adjusted_cns, theta2_adjusted_seg, theta2_subclonal_results, theta2_subclonal_cns,
      theta2_subclone_seg]
  
  expression_flatten_subclobal_results:
    run: ../tools/expression_flatten_list.cwl
    in:
      input_list: run_theta2_purity/theta2_subclonal_results
    out: [output]


$namespaces:
  sbg: https://sevenbridges.com
