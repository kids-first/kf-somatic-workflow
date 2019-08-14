cwlVersion: v1.0
class: Workflow
id: cnvkit-germline-wgs-wf

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  input_bam: { type: File, secondaryFiles: [.crai] }
  reference: {type: File, secondaryFiles: [.fai]}
  annotation_file: {type: File, doc: "refFlat.txt file"}
  output_basename: string

outputs:
  cnvkit_cnr: {type: File, outputSource: cnvkit/output_cnr}
  cnvkit_vcf: {type: File, outputSource: cnvkit/output_vcf}
  cnvkit_calls: {type: File, outputSource: cnvkit/output_calls}
  cnvkit_scatter: {type: File, outputSource: cnvkit/output_scatter}
  cnvkit_diagram: {type: File, outputSource: cnvkit/output_diagram}
  cnvkit_metrics: {type: File, outputSource: cnvkit/output_metrics}
  cnvkit_gainloss: {type: File, outputSource: cnvkit/output_gainloss}

steps:
  samtools_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_bam
      threads:
        valueFrom: ${return 36}
      reference: reference
    out: [bam_file]

  cnvkit: 
    run: ../tools/cnvkit_germline_wgs.cwl
    in:
      input_bam: samtools_cram2bam/bam_file
      reference: reference
      annotation_file: annotation_file
      output_basename: output_basename
    out: [output_cnr, output_vcf, output_calls, output_scatter, output_diagram, output_metrics, output_gainloss]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4
    