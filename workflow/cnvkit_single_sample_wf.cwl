cwlVersion: v1.0
class: Workflow
id: cnvkit-single-sample-wf

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  input_sample: { type: File, secondaryFiles: [.crai] }
  reference: {type: File, secondaryFiles: [.fai]}
  annotation_file: {type: File, doc: "refFlat.txt file"}
  output_basename: string
  wgs_mode: {type: ['null', string], doc: "for WGS mode, input Y. leave blank for hybrid mode"}
  capture_regions: {type: ['null', File], doc: "target regions for WES"}

outputs:
  cnvkit_cnr: {type: File, outputSource: cnvkit/output_cnr}
  cnvkit_vcf: {type: File, outputSource: cnvkit/output_vcf}
  cnvkit_calls: {type: File, outputSource: cnvkit/output_calls}
  cnvkit_scatter: {type: File, outputSource: cnvkit/output_scatter}
  cnvkit_diagram: {type: File, outputSource: cnvkit/output_diagram}
  cnvkit_metrics: {type: File, outputSource: cnvkit/output_metrics}
  cnvkit_gainloss: {type: File, outputSource: cnvkit/output_gainloss}

steps:
  samtools_sample_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_sample
      threads:
        valueFrom: ${return 36}
      reference: reference
    out: [bam_file]

  cnvkit: 
    run: ../tools/cnvkit_batch.cwl
    in:
      input_sample: samtools_sample_cram2bam/bam_file
      reference: reference
      annotation_file: annotation_file
      output_basename: output_basename
      wgs_mode: wgs_mode
      capture_regions: capture_regions
    out: [output_cnr, output_vcf, output_calls, output_scatter, output_diagram, output_metrics, output_gainloss]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4
    