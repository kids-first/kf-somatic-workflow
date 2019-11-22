cwlVersion: v1.0
class: Workflow
id: kfdrc_cnvkit_sub_wf

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  input_tumor_aligned: {type: File, secondaryFiles: ['^.bai']}
  input_normal_aligned: {type: File, secondaryFiles: ['^.bai']}
  reference: {type: File, secondaryFiles: [.fai]}
  b_allele_vcf: {type: ['null', File], doc: "b allele germline vcf, if available"}
  wgs_mode: {type: ['null', string], doc: "for WGS mode, input Y. leave blank for hybrid mode"}
  capture_regions: {type: ['null', File], doc: "target regions for WES"}
  annotation_file: {type: File, doc: "refFlat.txt file"}
  output_basename: string
  cnvkit_cnn_input: {type: ['null', File], doc: "If running using an existing .cnn, supply here"}
  threads: {type: ['null', int], default: 16}
  tumor_sample_name: {type: string, doc: "For seg file output and theta2 input"}
  normal_sample_name:  {type: string, doc: "For theta2 input"}
  sex: {type: string, doc: "Set sample sex.  CNVkit isn't always great at guessing it"}

outputs:
  cnvkit_cnr: {type: File, outputSource: cnvkit/output_cnr}
  cnvkit_cnn_output: {type: ['null', File], outputSource: cnvkit/output_cnn}
  cnvkit_calls: {type: File, outputSource: cnvkit/output_calls}
  cnvkit_metrics: {type: File, outputSource: cnvkit/output_metrics}
  cnvkit_gainloss: {type: File, outputSource: cnvkit/output_gainloss}
  cnvkit_seg: {type: File, outputSource: cnvkit/output_seg}

steps:

  cnvkit: 
    run: ../tools/cnvkit_batch.cwl
    in:
      input_sample: input_tumor_aligned
      input_control: input_normal_aligned
      reference: reference
      annotation_file: annotation_file
      output_basename: output_basename
      wgs_mode: wgs_mode
      capture_regions: capture_regions
      b_allele_vcf: b_allele_vcf
      threads: threads
      sex: sex
      tumor_sample_name: tumor_sample_name
      cnvkit_cnn: cnvkit_cnn_input
    out: [output_cnr, output_calls, output_scatter, output_diagram, output_metrics, output_gainloss, output_seg, output_cnn]
