cwlVersion: v1.0
class: Workflow
id: kfdrc_production_theta2_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
  - class: InlineJavascriptRequirement
inputs:
  # Required
  input_tumor_name: string
  input_normal_name: string
  min_theta2_frac: { type: 'float?', default: 0.01, doc: "Minimum fraction of genome with copy umber alterations.  Default is 0.05, recommend 0.01" }
  combined_include_expression: { type: 'string?', doc: "Theta2 Purity value: Filter expression if vcf has non-PASS combined calls, use as-needed, i.e. for VarDict: FILTER=\"PASS\" && (INFO/STATUS=\"Germline\" | INFO/STATUS=\"StrongSomatic\")" }
  combined_exclude_expression: { type: 'string?', doc: "Theta2 Purity value: Filter expression if vcf has non-PASS combined calls, use as-needed" }
  cnvkit_calls: { type: File, doc: "Calls from cnvkit" }
  cnvkit_cnn_output: { type: File, doc: "CNN calls from cnvkit" }
  vardict_prepass_vcf: { type: File, doc: "Prepass VCF output from Vardict" }

outputs:
  theta2_calls: { type: 'File?', outputSource: run_theta2_purity/theta2_adjusted_cns }
  theta2_seg: { type: 'File?', outputSource: run_theta2_purity/theta2_adjusted_seg }
  theta2_subclonal_results: { type: ['null', 'File[]'], outputSource: run_theta2_purity/theta2_subclonal_results }
  theta2_subclonal_cns: { type: ['null', 'File[]'], outputSource: run_theta2_purity/theta2_subclonal_cns }
  theta2_subclone_seg: { type: ['null', 'File[]'], outputSource: run_theta2_purity/theta2_subclone_seg }

steps:
  run_theta2_purity:
    run: ../sub_workflows/kfdrc_run_theta2_sub_wf.cwl
    in:
      tumor_cns: cnvkit_calls
      reference_cnn: cnvkit_cnn_output
      tumor_sample_name: input_tumor_name
      normal_sample_name: input_normal_name
      paired_vcf: vardict_prepass_vcf
      combined_include_expression: combined_include_expression
      combined_exclude_expression: combined_exclude_expression
      min_theta2_frac: min_theta2_frac
      output_basename:
        valueFrom: ${return inputs.tumor_cns.nameroot.split('.')[0]}
    out:
      [theta2_adjusted_cns, theta2_adjusted_seg, theta2_subclonal_results, theta2_subclonal_cns, theta2_subclone_seg]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 2
