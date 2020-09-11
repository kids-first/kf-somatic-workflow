cwlVersion: v1.0
class: Workflow
id: kfdrc_cnvkit_baf_theta2_wf

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  input_sample: {type: File}
  input_control: {type: ['null', File]}
  reference: {type: File, secondaryFiles: [.fai]}
  b_allele_vcf: {type: ['null', File], doc: "b allele germline vcf, if available"}
  paired_vcf: {type: File, doc: "Combined somatic and germline call file"}
  capture_regions: {type: ['null', File], doc: "target regions for WES"}
  annotation_file: {type: File, doc: "refFlat.txt file"}
  output_basename: string
  cnvkit_cnn_input: {type: ['null', File], doc: "If running using an existing .cnn, supply here"}
  wgs_mode: {type: ['null', string], doc: "for WGS mode, input Y. leave blank for hybrid mode", default: "Y"}
  threads: {type: ['null', int], default: 16}
  tumor_sample_name: {type: string, doc: "For seg file output and theta2 input"}
  normal_sample_name:  {type: string, doc: "For theta2 input"}
  sex: {type: string, doc: "Set sample sex.  CNVkit isn't always great at guessing it"}
  combined_include_expression: {type: ['null', string], doc: "Filter expression if vcf has non-PASS combined calls, use as-needed"}
  combined_exclude_expression: {type: ['null', string], doc: "Filter expression if vcf has non-PASS combined calls, use as-needed"}  
  min_theta2_frac: {type: ['null', float], doc: "Minimum fraction of genome with copy umber alterations.  Default is 0.05, recommend 0.01", default: 0.01}

outputs:
  cnvkit_cnr: {type: File, outputSource: cnvkit/output_cnr}
  cnvkit_cnn_output: {type: ['null', File], outputSource: cnvkit/output_cnn}
  cnvkit_calls: {type: File, outputSource: cnvkit/output_calls}
  cnvkit_metrics: {type: File, outputSource: cnvkit/output_metrics}
  cnvkit_gainloss: {type: File, outputSource: cnvkit/output_gainloss}
  cnvkit_seg: {type: File, outputSource: cnvkit/output_seg}
  theta2_calls: {type: File?, outputSource: cnvkit_import_theta2/theta2_adjusted_cns}
  theta2_seg: {type: File?, outputSource: cnvkit_import_theta2/theta2_adjusted_seg}
  theta2_subclonal_results: {type: ['null', 'File[]'], outputSource: [run_theta2/n3_graph, run_theta2/n2_results, run_theta2/best_results]}
  theta2_subclonal_cns: {type: ['null', 'File[]'], outputSource: cnvkit_import_theta2/theta2_subclone_cns}
  theta2_subclone_seg: {type: ['null', 'File[]'], outputSource: cnvkit_import_theta2/theta2_subclone_seg}

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
      tumor_cns: cnvkit/output_calls
      reference_cnn: cnvkit/output_cnn
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
      tumor_cns: cnvkit/output_calls
      tumor_sample_name: tumor_sample_name
      output_basename: output_basename
      theta2_best_results: run_theta2/best_results
      theta2_n2_results: run_theta2/n2_results
    out:
      [theta2_adjusted_cns, theta2_adjusted_seg, theta2_subclone_cns, theta2_subclone_seg]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 2
    