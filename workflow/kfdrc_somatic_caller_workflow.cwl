cwlVersion: v1.0
class: Workflow
id: kfdrc_somatic_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict]}
  reference_dict: File
  wgs_calling_interval_list: File
  hg38_strelka_bed: File
  input_tumor_aligned: File
  input_normal_aligned: File
  exome_flag: string
  vep_cache: {type: File, label: tar gzipped cache from ensembl/local converted cache}
  output_basename: string

outputs:
  vardict_vep_vcf: {type: File, outputSource: vep_annot_vardict/output_vcf}
  vardict_vep_html: {type: File, outputSource: vep_annot_vardict/output_html}
  strelka2_vep_vcf: {type: File, outputSource: vep_annot_strelka2/output_vcf}
  strelka2_vep_html: {type: File, outputSource: vep_annot_strelka2/output_html}

steps:
  gatk_intervallisttools:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: wgs_calling_interval_list
    out: [output]

  strelka2:
    in:
      input_tumor_aligned: input_tumor_aligned
      input_normal_aligned: input_normal_aligned
      reference: indexed_reference_fasta
      hg38_strelka_bed: hg38_strelka_bed
      exome_flag: exome_flag
    out: [output]
    run: ../tools/strelka2.cwl

  vardict_somatic:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.4xlarge;ebs-gp2;500
    run: ../tools/vardictjava.cwl
    in:
      input_tumor_bam: input_tumor_aligned
      input_normal_bam: input_normal_aligned
      bed: gatk_intervallisttools/output
      reference: indexed_reference_fasta
      output_basename: output_basename
    scatter: [interval_list]
    out: [vardict_vcf]
  
  merge_vardict_vcf:
    run: ../tools/gatk_mergevcfs_pass_filter.cwl
    label: Merge & pass filter VarDict
    in:
      input_vcfs: vardict_somatic/output
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name:
        valueFrom: ${ return "vardict"}
    out: [merged_vcf]

  merge_strelka2_vcf:
    run: ../tools/gatk_mergevcfs_pass_filter.cwl
    label: Merge & pass filter strekla2
    in:
      input_vcfs: [strelka2/output_snv, strelka2/output_indel]
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name:
        valueFrom: ${ return "strelka2"}
    out: [merged_vcf]
  
  vep_annot_vardict:
    run: ../tools/variant_effect_predictor.cwl
    in:
      input_vcf: merge_vardict_vcf/merged_vcf
      output_basename: output_basename
      tool_name:
        valueFrom: ${ return "vardict"}
      reference: indexed_reference_fasta
      cache: vep_cache
    out: [output_vcf, output_html, warn_txt]
  
  vep_annot_strelka2:
    run: ../tools/variant_effect_predictor.cwl
    in:
      input_vcf: merge_strelka2_vcf/merged_vcf
      output_basename: output_basename
      tool_name:
        valueFrom: ${ return "strelka2"}
      reference: indexed_reference_fasta
      cache: vep_cache
    out: [output_vcf, output_html, warn_txt]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 3