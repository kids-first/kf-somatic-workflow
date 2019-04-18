cwlVersion: v1.0
class: Workflow
id: vardict_test_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict]}
  wgs_calling_interval_list: File
  input_tumor_aligned: File
  input_tumor_name: string
  input_normal_aligned: File
  input_normal_name: string
  output_basename: string
  reference_dict: File

outputs:
  vardict_vcf: {type: File, outputSource: merge_vardict_vcf/merged_vcf}

steps:
  gatk_intervallisttools:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: wgs_calling_interval_list
    out: [output]

  vardict:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: m5.4xlarge
    run: ../tools/vardictjava.cwl
    in:
      input_tumor_bam: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_bam: input_normal_aligned
      input_normal_name: input_normal_name
      reference: indexed_reference_fasta
      bed: gatk_intervallisttools/output
      output_basename: output_basename
    scatter: [bed]
    out: [vardict_vcf]

  merge_vardict_vcf:
    run: ../tools/gatk_mergevcfs.cwl
    label: GATK Merge vardict
    in:
      input_vcfs: vardict/vardict_vcf
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name:
        valueFrom: ${return "vardict"}
      silent_flag:
        valueFrom: ${return 1}
    out: [merged_vcf]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4