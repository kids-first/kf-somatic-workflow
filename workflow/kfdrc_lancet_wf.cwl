cwlVersion: v1.0
class: Workflow
id: kfdrc_lancet_wf
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
  exome_flag: {type: ['null', string], doc: "set to 'Y' for exome mode"}
  select_vars_mode: {type: string, doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression"}
  vep_cache: {type: File, label: tar gzipped cache from ensembl/local converted cache}
  window: { type: int, doc: "window size for lancet.  default is 600"}

outputs:
  lancet_vep_vcf: {type: File, outputSource: vep_annot_lancet/output_vcf}
  lancet_vep_tbi: {type: File, outputSource: vep_annot_lancet/output_tbi}
  lancet_vep_maf: {type: File, outputSource: vep_annot_lancet/output_maf}
  lancet_prepass_vcf: {type: File, outputSource: sort_merge_lancet_vcf/merged_vcf}

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
    run: ../tools/lancet.cwl
    in:
      input_tumor_bam: samtools_tumor_cram2bam/bam_file
      input_normal_bam: samtools_normal_cram2bam/bam_file
      reference: indexed_reference_fasta
      bed: gatk_intervallisttools/output
      output_basename: output_basename
      window: window
    scatter: [bed]
    out: [lancet_vcf]

  sort_merge_lancet_vcf:
    run: ../dev/gatk_sortvcf.cwl
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
    run: ../tools/gatk_selectvariants.cwl
    label: GATK Select Lancet PASS
    in:
      input_vcf: sort_merge_lancet_vcf/merged_vcf
      output_basename: output_basename
      tool_name:
        valueFrom: ${return "lancet"}
      mode: select_vars_mode
    out: [pass_vcf]

  vep_annot_lancet:
    run: ../tools/vep_vcf2maf.cwl
    in:
      input_vcf: gatk_selectvariants_lancet/pass_vcf
      output_basename: output_basename
      tumor_id: input_tumor_name
      normal_id: input_normal_name
      tool_name:
        valueFrom: ${return "lancet_somatic"}
      reference: indexed_reference_fasta
      cache: vep_cache
    out: [output_vcf, output_tbi, output_maf, warn_txt]


$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4