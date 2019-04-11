cwlVersion: v1.0
class: Workflow
id: mutect2_support
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict]}
  reference_dict: File
  wgs_calling_interval_list: File
  input_tumor_aligned: File
  input_normal_aligned: File
  exac_common_vcf: {type: File, secondaryFiles: [.tbi]}
  output_basename: string
  f1r2_counts: {type: "File[]", doc: "orientation counts from mutect2 outputs"}

outputs:
  merged_pileup_summary: {type: File, outputSource: gather_pileup_summaries/merged_table}
  
steps:

  gatk_learn_orientation_bias:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    run: ../tools/gatk_mergevcfs_pass_filter.cwl
    label: Merge mutect2 output
    in:
      input_tgz: f1r2_counts
      output_basename: output_basename
      tool_name:
        valueFrom: ${return "mutect2"}
    out: [f1r2_bias]

  gatk_intervallisttools:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: wgs_calling_interval_list
    out: [output]

  gatk_get_tumor_pileup_summaries:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge;ebs-gp2;500
    run: ../tools/gatk_getpileupsummaries.cwl
    in:
      aligned_reads: input_tumor_aligned
      reference: indexed_reference_fasta
      interval_list: gatk_intervallisttools/output
      exac_common_vcf: exac_common_vcf
    scatter: [interval_list]
    out: [pileup_table]

  gatk_get_normal_pileup_summaries:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge;ebs-gp2;500
    run: ../tools/gatk_getpileupsummaries.cwl
    in:
      aligned_reads: input_normal_aligned
      reference: indexed_reference_fasta
      interval_list: gatk_intervallisttools/output
      exac_common_vcf: exac_common_vcf
    scatter: [interval_list]
    out: [pileup_table]

  gatk_gather_tumor_pileup_summaries:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    run: ../tools/gatk_gatherpileupsummaries.cwl
    label: GATK merge pileup tables
    in:
      input_tables: gatk_get_tumor_pileup_summaries/pileup_table
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name:
        valueFrom: ${return "mutect2"}
    out: [merged_table]

  gatk_gather_normal_pileup_summaries:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    run: ../tools/gatk_gatherpileupsummaries.cwl
    label: GATK merge pileup tables
    in:
      input_tables: gatk_get_normal_pileup_summaries/pileup_table
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name:
        valueFrom: ${return "mutect2"}
    out: [merged_table]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 2