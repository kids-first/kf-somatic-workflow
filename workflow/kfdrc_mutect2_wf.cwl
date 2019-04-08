cwlVersion: v1.0
class: Workflow
id: kfdrc_mutect2_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict]}
  reference_dict: File
  wgs_calling_interval_list: File
  input_tumor_aligned: File
  input_tumor_name: string
  input_normal_aligned: File
  input_normal_name: string
  exome_flag: ['null', string]
  vep_cache: {type: File, label: tar gzipped cache from ensembl/local converted cache}
  output_basename: string

outputs:
  mutect2_vep_vcf: {type: File, outputSource: vep_annot_mutect2/output_vcf}
  mutect2_vep_tbi: {type: File, outputSource: vep_annot_mutect2/output_tbi}
  mutect2_vep_maf: {type: File, outputSource: vep_annot_mutect2/output_maf}
  
steps:
  gatk_intervallisttools:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: wgs_calling_interval_list
    out: [output]

  mutect2:
    run: ../tools/gatk_Mutect2.cwl
    in:
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      reference: indexed_reference_fasta
      interval_list: gatk_intervallisttools/output
      exome_flag: exome_flag
    scatter: [interval_list]
    out: [mutect2_vcf]

  merge_mutect2_vcf:
    run: ../tools/gatk_mergevcfs_pass_filter.cwl
    label: Merge & pass filter mutect2
    in:
      input_vcfs: mutect2/mutect2_vcf
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name:
        valueFrom: ${return "mutect2"}
    out: [merged_vcf]

  vep_annot_mutect2:
    run: ../tools/vep_vcf2maf.cwl
    in:
      input_vcf: merge_mutect2_vcf/merged_vcf
      output_basename: output_basename
      tumor_id: input_tumor_name
      normal_id: input_normal_name
      tool_name:
        valueFrom: ${return "mutect2_somatic"}
      reference: indexed_reference_fasta
      cache: vep_cache
    out: [output_vcf, output_tbi, output_html, warn_txt]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4
  - class: 'sbg:AWSInstanceType'
    value: c5.9xlarge;ebs-gp2;750
