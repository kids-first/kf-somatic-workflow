cwlVersion: v1.0
class: Workflow
id: kfdrc_cnvkit_batch_wf

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  input_sample: { type: File, secondaryFiles: [.crai] }
  input_control: { type: ['null', File], secondaryFiles: [.crai] }
  reference: {type: File, secondaryFiles: [.fai]}
  b_allele_vcf: {type: ['null', File], doc: "b allele germline vcf, if available"}
  capture_regions: {type: ['null', File], doc: "target regions for WES"}
  annotation_file: {type: File, doc: "refFlat.txt file"}
  output_basename: string
  cnvkit_cnn_input: {type: ['null', File], doc: "If running using an existing .cnn, supply here"}
  wgs_mode: {type: ['null', string], doc: "for WGS mode, input Y. leave blank for hybrid mode"}
  threads: {type: ['null', int], default: 16}
  tumor_sample_name: {type: string, doc: "For seg file output"}
  sex: {type: string, doc: "Set sample sex.  CNVkit isn't always great at guessing it"}

outputs:
  cnvkit_cnr: {type: File, outputSource: cnvkit/output_cnr}
  cnvkit_cnn_output: {type: ['null', File], outputSource: cnvkit/output_cnn}
  cnvkit_calls: {type: File, outputSource: cnvkit/output_calls}
  cnvkit_scatter: {type: File, outputSource: cnvkit/output_scatter}
  cnvkit_diagram: {type: File, outputSource: cnvkit/output_diagram}
  cnvkit_metrics: {type: File, outputSource: cnvkit/output_metrics}
  cnvkit_gainloss: {type: File, outputSource: cnvkit/output_gainloss}
  cnvkit_seg: {type: File, outputSource: cnvkit/output_seg}

steps:
  gatk_filter_germline:
    run: ../tools/gatk_filter_germline_variant.cwl
    in:
      input_vcf: b_allele_vcf
      reference_fasta: reference
      output_basename: output_basename
    out:
      [filtered_vcf, filtered_pass_vcf]

  samtools_sample_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_sample
      reference: reference
    out: [bam_file]

  samtools_control_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_control
      reference: reference
    out: [bam_file]

  cnvkit: 
    run: ../tools/cnvkit_batch.cwl
    in:
      input_sample: samtools_sample_cram2bam/bam_file
      input_control: samtools_control_cram2bam/bam_file
      reference: reference
      annotation_file: annotation_file
      output_basename: output_basename
      wgs_mode: wgs_mode
      capture_regions: capture_regions
      b_allele_vcf: gatk_filter_germline/filtered_pass_vcf
      threads: threads
      sex: sex
      tumor_sample_name: tumor_sample_name
      cnvkit_cnn: cnvkit_cnn_input
    out: [output_cnr, output_calls, output_scatter, output_diagram, output_metrics, output_gainloss, output_seg, output_cnn]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 2
    