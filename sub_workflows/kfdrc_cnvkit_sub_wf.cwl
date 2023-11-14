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
  input_normal_aligned: {type: File, secondaryFiles: ['^.bai']}
  reference: {type: File, secondaryFiles: [.fai]}
  b_allele_vcf: {type: ['null', File], doc: "b allele germline vcf, if available"}
  wgs_mode: {type: ['null', string], doc: "for WGS mode, input Y. leave blank for hybrid mode"}
  capture_regions: {type: ['null', File], doc: "target regions for WES. These are the bait regions."}
  blacklist_regions: {type: 'File?', doc: "Regions of the genome over which CNVkit should not perform resequencing. CNVkit refers to these regions as inaccessible." }
  annotation_file: {type: File, doc: "refFlat.txt file"}
  output_basename: string
  cnvkit_cnn_input: {type: ['null', File], doc: "If running using an existing .cnn, supply here"}
  threads: {type: ['null', int], default: 16}
  tumor_sample_name: {type: string, doc: "For seg file output and theta2 input"}
  normal_sample_name:  {type: string, doc: "For theta2 input"}
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
    run: ../tools/cnvkit_access.cwl
    in:
      reference_fasta: reference
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
      wgs_mode: wgs_mode
      access: access/bed
      capture_regions: capture_regions
      b_allele_vcf: b_allele_vcf
      threads: threads
      sex: sex
      tumor_sample_name: tumor_sample_name
      cnvkit_cnn: cnvkit_cnn_input
    out: [output_cnr, output_calls, output_cns, output_scatter, output_diagram, output_metrics, output_gainloss, output_seg, output_cnn]
