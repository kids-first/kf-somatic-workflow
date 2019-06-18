cwlVersion: v1.0
class: Workflow
id: kfdrc_lancet_wf_benchmark
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict]}
  wgs_calling_interval_list: File
  input_tumor_aligned: {type: File, secondaryFiles: [^.bai]}
  input_normal_aligned: {type: File, secondaryFiles: [^.bai]}
  output_basename: string
  reference_dict: File
  exome_flag: {type: ['null', string], doc: "set to 'Y' for exome mode"}
  select_vars_mode: {type: string, doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression"}
  window: { type: int, doc: "window size for lancet.  default is 600"}

outputs:
  lancet_pass_vcf: {type: File, outputSource: gatk_selectvariants_lancet/pass_vcf}
  lancet_prepass_vcf: {type: File, outputSource: sort_merge_lancet_vcf/merged_vcf}

steps:
  gatk_intervallisttools:
    run: ../../tools/gatk_intervallisttool.cwl
    in:
      interval_list: wgs_calling_interval_list
      reference_dict: reference_dict
      exome_flag: exome_flag
      scatter_ct:
        valueFrom: ${return 50}
      bands:
        valueFrom: ${return 80000000}
    out: [output]

  lancet:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge;ebs-gp2;500
    run: ../../tools/lancet.cwl
    in:
      input_tumor_bam: input_tumor_aligned
      input_normal_bam: input_normal_aligned
      reference: indexed_reference_fasta
      bed: gatk_intervallisttools/output
      output_basename: output_basename
      window: window
    scatter: [bed]
    out: [lancet_vcf]

  sort_merge_lancet_vcf:
    run: ../../tools/gatk_sortvcf.cwl
    label: GATK Sort & Merge lancet
    in:
      input_vcfs: lancet/lancet_vcf
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name:
        valueFrom: ${return "lancet"}
    out: [merged_vcf]

  gatk_selectvariants_lancet:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    run: ../../tools/gatk_selectvariants.cwl
    label: GATK Select Lancet PASS
    in:
      input_vcf: sort_merge_lancet_vcf/merged_vcf
      output_basename: output_basename
      tool_name:
        valueFrom: ${return "lancet"}
      mode: select_vars_mode
    out: [pass_vcf]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4