cwlVersion: v1.0
class: Workflow
id: temp_rerun_vcf2maf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  reference: {type: File,  secondaryFiles: [.fai], label: Fasta genome assembly with index}
  input_vcf:
    type: File
    secondaryFiles: [.tbi]
  output_basename: string
  tumor_id: string
  normal_id: string
  tool_name: string
  cache: {type: File, label: tar gzipped cache from ensembl/local converted cache}
  strip_info: {type: ['null', string], doc: "If given, remove previous annotation information based on INFO file, i.e. to strip VEP info, use INFO/ANN"}

outputs:
  annotated_vcf: {type: File, outputSource: vep_annot_vcf/output_vcf}
  annotated_tbi: {type: File, outputSource: vep_annot_vcf/output_tbi}
  annotated_maf: {type: File, outputSource: vep_annot_vcf/output_maf}

steps:
  bcftools_strip_ann:
    run: ../dev/bcftools_strip_ann.cwl
    in:
      input_vcf: input_vcf
      output_basename: output_basename
      tool_name: tool_name
      strip_info: strip_info
    out:
      [stripped_vcf]

  vep_annot_vcf:
    run: ../tools/vep_vcf2maf.cwl
    in:
      input_vcf: bcftools_strip_ann/stripped_vcf
      output_basename: output_basename
      tumor_id: tumor_id
      normal_id: normal_id
      tool_name: tool_name
      reference: reference
      cache: cache
    out: [output_vcf, output_tbi, output_maf, warn_txt]