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
  input_tumor_name: string
  input_normal_aligned: File
  input_normal_name: string
  exome_flag: ['null', string]
  vep_cache: {type: File, label: tar gzipped cache from ensembl/local converted cache}
  output_basename: string

outputs:
  vardict_vep_vcf: {type: File, outputSource: vep_annot_vardict/output_vcf}
  vardict_vep_html: {type: File, outputSource: vep_annot_vardict/output_html}
  strelka2_vep_vcf: {type: File, outputSource: vep_annot_strelka2/output_vcf}
  strelka2_vep_html: {type: File, outputSource: vep_annot_strelka2/output_html}
  mutect2_vep_vcf: {type: File, outputSource: vep_annot_mutect2/output_vcf}
  mutect2_vep_html: {type: File, outputSource: vep_annot_mutect2/output_html}

steps:
  samtools_tumor_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_tumor_aligned
      threads:
        valueFrom: ${return 36}
      reference: indexed_reference_fasta
    out: [bam_file]

  samtools_normal_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_normal_aligned
      threads:
        valueFrom: ${return 36}
      reference: indexed_reference_fasta
    out: [bam_file]

  gatk_intervallisttools:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: wgs_calling_interval_list
    out: [output]

  strelka2:
    run: ../tools/strelka2.cwl
    in:
      input_tumor_aligned: input_tumor_aligned
      input_normal_aligned: input_normal_aligned
      reference: indexed_reference_fasta
      hg38_strelka_bed: hg38_strelka_bed
      exome_flag: exome_flag
    out: [output]

  vardict_somatic:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.4xlarge;ebs-gp2;500
    run: ../tools/vardictjava.cwl
    in:
      input_tumor_bam: samtools_tumor_cram2bam/bam_file
      input_tumor_name: input_tumor_name
      input_normal_bam: samtools_normal_cram2bam/bam_file
      input_normal_name: input_normal_name
      bed: gatk_intervallisttools/output
      reference: indexed_reference_fasta
      output_basename: output_basename
    scatter: [bed]
    out: [vardict_vcf]

  mutect2:
    run: ../tools/gatk_Mutect2.cwl
    in:
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      reference: indexed_reference_fasta
      interval_list: gatk_intervallisttools/output
    scatter: [interval_list]
    out: [mutect2_vcf]
  
  merge_vardict_vcf:
    run: ../tools/gatk_mergevcfs_pass_filter.cwl
    label: Merge & pass filter VarDict
    in:
      input_vcfs: vardict_somatic/vardict_vcf
      input_tumor_name: input_tumor_name
      input_normal_name: input_normal_name
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

  rename_strelka_samples:
    run: ../tools/bcftools_reheader_vcf.cwl
    in:
      input_vcf: merge_strelka2_vcf/merged_vcf
      input_normal_name: input_normal_name
      input_tumor_name: input_tumor_name
    out: [reheadered_vcf]

  merge_mutect2_vcf:
    run: ../tools/gatk_mergevcfs_pass_filter.cwl
    label: Merge & pass filter mutect2
    in:
      input_vcfs: mutect2/mutect2_vcf
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name:
        valueFrom: ${ return "mutect"}
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
      input_vcf: rename_strelka_samples/reheadered_vcf
      output_basename: output_basename
      tool_name:
        valueFrom: ${ return "strelka2"}
      reference: indexed_reference_fasta
      cache: vep_cache
    out: [output_vcf, output_html, warn_txt]

  vep_annot_mutect2:
    run: ../tools/variant_effect_predictor.cwl
    in:
      input_vcf: merge_mutect2_vcf/merged_vcf
      output_basename: output_basename
      tool_name:
        valueFrom: ${ return "mutect2"}
      reference: indexed_reference_fasta
      cache: vep_cache
    out: [output_vcf, output_html, warn_txt]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 3