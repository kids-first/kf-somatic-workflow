cwlVersion: v1.2
class: Workflow
id: kfdrc_cnvkit_sub_wf

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:
  input_tumor_aligned: {type: File, secondaryFiles: ['^.bai']}
  input_normal_aligned: {type: 'File?', secondaryFiles: ['^.bai'],
    doc: "Only needed if cnn file not available"}
  reference: {type: File, secondaryFiles: [.fai]}
  b_allele_vcf: {type: ['null', File], doc: "b allele germline vcf, if available"}
  run_mode: { type: ['null', {type: enum, name: run_mode, symbols: ["hybrid", "wgs", "amplicon"]}],
    default: "hybrid",
    doc: "for WGS mode, input wgs. leave blank for hybrid mode" }
  target_regions: {type: ['null', File], doc: "target regions for analysis. required for WES, useful if already masked for WGS"}
  blacklist_regions: {type: 'File?',
    doc: "Regions of the genome over which CNVkit should not perform resequencing. CNVkit refers to these regions as inaccessible. Skipped if target_regions given" }
  annotation_file: {type: 'File?', doc: "refFlat.txt file. Needed if cnn doesn't already exist"}
  output_basename: string
  cnvkit_cnn_input: {type: ['null', File], doc: "If running using an existing .cnn, supply here"}
  threads: {type: ['null', int], default: 16}
  tumor_sample_name: {type: string, doc: "For seg file output and theta2 input"}
  sex: {type: ['null', {type: enum, name: sex, symbols: ["x", "y"]}], doc: "Sex, for simplicity x for female y for male", default: "x"}

outputs:
  cnvkit_cnr: {type: File, outputSource: batch/output_cnr}
  cnvkit_cnn_output: {type: ['null', File], outputSource: batch/output_cnn}
  cnvkit_cns: { type: File, outputSource: batch/output_cns }
  cnvkit_calls: {type: File, outputSource: batch/output_calls}
  cnvkit_metrics: {type: File, outputSource: batch/output_metrics}
  cnvkit_gainloss: {type: File, outputSource: batch/output_gainloss}
  cnvkit_seg: {type: File, outputSource: batch/output_seg}
  cnvkit_scatter_plot: {type: File, outputSource: batch/output_scatter}
  cnvkit_diagram: {type: File, outputSource: batch/output_diagram}

steps:
  access:
    when: $(inputs.target_regions == null && inputs.blacklist_regions)
    run: ../tools/cnvkit_access.cwl
    in:
      reference_fasta: reference
      target_regions: target_regions
      min_gap_size:
        valueFrom: $(200)
      exclude: blacklist_regions
    out: [bed]
  batch:
    run: ../tools/cnvkit_batch.cwl
    in:
      input_sample: input_tumor_aligned
      input_control: input_normal_aligned
      reference: reference
      annotation_file: annotation_file
      output_basename: output_basename
      run_mode: run_mode
      access: access/bed
      target_regions:
        source: [cnvkit_cnn_input, target_regions]
        valueFrom: |
          $(self[0] == null ? self[1] : null)
      b_allele_vcf: b_allele_vcf
      threads: threads
      sex: sex
      tumor_sample_name: tumor_sample_name
      cnvkit_cnn: cnvkit_cnn_input
    out: [output_cnr, output_calls, output_cns, output_scatter, output_diagram, output_metrics, output_gainloss, output_seg, output_cnn]
