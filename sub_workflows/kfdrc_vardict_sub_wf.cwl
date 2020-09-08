cwlVersion: v1.0
class: Workflow
id: kfdrc_vardict_1_7_sub_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict]}
  input_tumor_aligned: {type: File, secondaryFiles: ['^.bai']}
  input_tumor_name: string
  input_normal_aligned: {type: File, secondaryFiles: ['^.bai']}
  input_normal_name: string
  output_basename: string
  reference_dict: File
  padding: {type: ['null', int], doc: "Padding to add to input intervals, recommened 0 if intervals already padded, 150 if not", default: 150}
  bed_invtl_split: {type: 'File[]', doc: "Bed file intervals passed on from and outside pre-processing step"}
  min_vaf: {type: ['null', float], doc: "Min variant allele frequency for vardict to consider.  Recommend 0.05", default: 0.05}
  select_vars_mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression", default: "gatk"}
  cpus: {type: ['null', int], default: 9}
  ram: {type: ['null', int], default: 18, doc: "In GB"}
  vep_cache: {type: File, label: tar gzipped cache from ensembl/local converted cache}
  vep_ref_build: {type: ['null', string], doc: "Genome ref build used, should line up with cache.", default: "GRCh38" }

outputs:
  vardict_vep_somatic_only_vcf: {type: File, outputSource: vep_annot_vardict/output_vcf}
  vardict_vep_somatic_only_tbi: {type: File, outputSource: vep_annot_vardict/output_tbi}
  vardict_vep_somatic_only_maf: {type: File, outputSource: vep_annot_vardict/output_maf}
  vardict_prepass_vcf: {type: File, outputSource: sort_merge_vardict_vcf/merged_vcf}

steps:

  vardict:
    run: ../tools/vardictjava.cwl
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge
    in:
      input_tumor_bam: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_bam: input_normal_aligned
      input_normal_name: input_normal_name
      padding: padding
      min_vaf: min_vaf
      cpus: cpus
      ram: ram
      reference: indexed_reference_fasta
      bed: bed_invtl_split
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
    run: ../tools/bcbio_filter_vardict_somatic.cwl
    in:
      input_vcf: sort_merge_vardict_vcf/merged_vcf
      output_basename: output_basename
    out: [filtered_vcf] 

  gatk_selectvariants_vardict:
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
      ref_build: vep_ref_build
    out: [output_vcf, output_tbi, output_maf, warn_txt]

$namespaces:
  sbg: https://sevenbridges.com
