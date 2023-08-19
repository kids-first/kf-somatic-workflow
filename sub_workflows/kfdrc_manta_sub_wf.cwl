cwlVersion: v1.2
class: Workflow
id: kfdrc_manta_sub_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  indexed_reference_fasta: {type: 'File', secondaryFiles: [{ pattern: '.fai', required: true }, { pattern: '^.dict', required: true }]}
  reference_dict: File
  hg38_strelka_bed: {type: 'File', secondaryFiles: [{ pattern: '.tbi', required: true }]}
  input_tumor_aligned:
    type: File
    secondaryFiles: [{ pattern: '.bai', required: false }, { pattern: '^.bai', required: false }, { pattern: '.crai', required: false }, { pattern: '^.crai', required: false }]
    doc: "tumor BAM or CRAM"
  input_tumor_name: string
  old_tumor_name: {type: 'string?', doc: "If `SM:` sample name in te align file is different than `input_tumor_name`, you **must** provide it here"}
  input_normal_aligned:
    type: File
    secondaryFiles: [{ pattern: '.bai', required: false }, { pattern: '^.bai', required: false }, { pattern: '.crai', required: false }, { pattern: '^.crai', required: false }]
    doc: "normal BAM or CRAM"
  input_normal_name: string
  old_normal_name: {type: 'string?', doc: "If `SM:` sample name in te align file is different than `input_normal_name`, you **must** provide it here"}
  vep_cache: {type: 'File', label: tar gzipped cache from ensembl/local converted cache}
  output_basename: string
  manta_memory: {type: 'int?'}
  manta_cores: {type: 'int?'}
  select_vars_mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression", default: "gatk"}

outputs:
  manta_prepass_vcf: {type: 'File', outputSource: pickvalue_workaround/output}
  manta_pass_vcf: {type: 'File', outputSource: gatk_selectvariants_manta/pass_vcf}
  manta_small_indels: {type: 'File', outputSource: manta/small_indels}

steps:
  manta:
    run: ../tools/manta.cwl
    in:
      input_tumor_cram: input_tumor_aligned
      input_normal_cram: input_normal_aligned
      output_basename: output_basename
      ram: manta_memory
      cores: manta_cores
      reference: indexed_reference_fasta
      hg38_strelka_bed: hg38_strelka_bed
    out: [output_sv, small_indels]

  rename_manta_samples:
    run: ../tools/bcftools_reheader_samples_index.cwl
    when: $(inputs.old_tumor_name != null && inputs.old_tumor_name != null)
    in:
      input_vcf: manta/output_sv
      output_filename:
        valueFrom: |
          $(inputs.input_vcf.basename.replace(".vcf", ".reheadered.vcf"))
      new_normal_name: input_normal_name
      new_tumor_name: input_tumor_name
      old_normal_name: old_normal_name
      old_tumor_name: old_tumor_name
      tbi:
        valueFrom: |
          $(1 == 1)
    out: [reheadered_vcf]

  pickvalue_workaround:
    run: ../tools/expression_pickvalue_workaround.cwl
    in:
      input_file:
        source: [rename_manta_samples/reheadered_vcf, manta/output_sv]
        pickValue: first_non_null
    out: [output]

  gatk_selectvariants_manta:
    run: ../tools/gatk_selectvariants.cwl
    label: GATK Select Manta PASS
    in:
      input_vcf: pickvalue_workaround/output
      output_basename: output_basename
      tool_name:
        valueFrom: ${return "manta"}
      mode: select_vars_mode
    out: [pass_vcf]

