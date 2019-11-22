cwlVersion: v1.0
class: Workflow
id: kfdrc_vardict_1_7_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict]}
  wgs_calling_interval_list: {type: File, doc: "WGS bed file with N's regions removed"}
  input_tumor_aligned: File
  input_tumor_name: string
  input_normal_aligned: File
  input_normal_name: string
  output_basename: string
  reference_dict: File
  min_vaf: {type: ['null', float], doc: "Min variant allele frequency for vardict to consider.  Recommend 0.05", default: 0.05}
  select_vars_mode: {type: string, doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression"}
  cpus: {type: ['null', int], default: 9}
  ram: {type: ['null', int], default: 18, doc: "In GB"}
  vep_cache: {type: File, label: tar gzipped cache from ensembl/local converted cache}

outputs:
  vardict_vep_somatic_only_vcf: {type: File, outputSource: vep_annot_vardict/output_vcf}
  vardict_vep_somatic_only_tbi: {type: File, outputSource: vep_annot_vardict/output_tbi}
  vardict_vep_somatic_only_maf: {type: File, outputSource: vep_annot_vardict/output_maf}
  vardict_prepass_vcf: {type: File, outputSource: sort_merge_vardict_vcf/merged_vcf}

steps:
  samtools_tumor_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_tumor_aligned
      threads:
        valueFrom: ${return 16}
      reference: indexed_reference_fasta
    out: [bam_file]
  samtools_normal_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_normal_aligned
      threads:
        valueFrom: ${return 16}
      reference: indexed_reference_fasta
    out: [bam_file]

  python_vardict_interval_split:
    run: ../tools/python_vardict_interval_split.cwl
    in:
      wgs_bed_file: wgs_calling_interval_list
    out: [split_intervals_bed]

  vardict:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge;ebs-gp2;750
    run: ../tools/vardictjava.cwl
    in:
      input_tumor_bam: samtools_tumor_cram2bam/bam_file
      input_tumor_name: input_tumor_name
      input_normal_bam: samtools_normal_cram2bam/bam_file
      input_normal_name: input_normal_name
      min_vaf: min_vaf
      cpus: cpus
      ram: ram
      reference: indexed_reference_fasta
      bed: python_vardict_interval_split/split_intervals_bed
      output_basename: output_basename
    scatter: [bed]
    out: [vardict_vcf]

  sort_merge_vardict_vcf:
    run: ../tools/gatk_sortvcf.cwl
    label: GATK Sort & merge vardict
    in:
      input_vcfs: vardict/vardict_vcf
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name:
        valueFrom: ${return "vardict"}
    out: [merged_vcf]

  bcbio_filter_fp_somatic:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    run: ../tools/bcbio_filter_vardict_somatic.cwl
    in:
      input_vcf: sort_merge_vardict_vcf/merged_vcf
      output_basename: output_basename
    out: [filtered_vcf] 

  gatk_selectvariants_vardict:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    run: ../tools/gatk_selectvariants.cwl
    label: GATK Select Vardict PASS
    in:
      input_vcf: bcbio_filter_fp_somatic/filtered_vcf
      output_basename: output_basename
      tool_name:
        valueFrom: ${return "vardict"}
      mode: select_vars_mode
    out: [pass_vcf]

  vep_annot_vardict:
    run: ../tools/vep_vcf2maf.cwl
    in:
      input_vcf: gatk_selectvariants_vardict/pass_vcf
      output_basename: output_basename
      tumor_id: input_tumor_name
      normal_id: input_normal_name
      tool_name:
        valueFrom: ${return "vardict_somatic"}
      reference: indexed_reference_fasta
      cache: vep_cache
    out: [output_vcf, output_tbi, output_maf, warn_txt]


$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4