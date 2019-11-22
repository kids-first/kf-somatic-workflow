cwlVersion: v1.0
class: Workflow
id: kfdrc_consensus_calling
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai]}
  strelka2_vcf: {type: File, secondaryFiles: ['.tbi']}
  mutect2_vcf: {type: File, secondaryFiles: ['.tbi']}
  lancet_vcf: {type: File, secondaryFiles: ['.tbi']}
  vardict_vcf: {type: File, secondaryFiles: ['.tbi']}
  output_basename: string
  min_overlap: {type: ['null', int], default: 2, doc: "Min number of callers to declare consensus.  Default is 2"}
  strip_info: {type: ['null', string], doc: "If given, remove previous annotation information based on INFO file, i.e. to strip VEP info, use INFO/ANN"}

outputs:
  consensus_vcf: {type: File, outputSource: bcbio_variant_recall_ensemble/consensus_vcf}

steps:
  normalize_strelka2_vcf:
    run: ../tools/normalize_vcf.cwl
    in:
      input_vcf: strelka2_vcf
      indexed_reference_fasta: indexed_reference_fasta
      output_basename: output_basename
      tool_name:
        valueFrom: ${return "strelka2";}
      strip_info: strip_info
    out: [normalized_vcf]

  normalize_mutect2_vcf:
    run: ../tools/normalize_vcf.cwl
    in:
      input_vcf: mutect2_vcf
      indexed_reference_fasta: indexed_reference_fasta
      output_basename: output_basename
      tool_name:
        valueFrom: ${return "mutect2";}
      strip_info: strip_info
    out: [normalized_vcf]

  normalize_lancet_vcf:
    run: ../tools/normalize_vcf.cwl
    in:
      input_vcf: lancet_vcf
      indexed_reference_fasta: indexed_reference_fasta
      output_basename: output_basename
      tool_name:
        valueFrom: ${return "lancet";}
      strip_info: strip_info
    out: [normalized_vcf]

  normalize_vardict_vcf:
    run: ../tools/normalize_vcf.cwl
    in:
      input_vcf: vardict_vcf
      indexed_reference_fasta: indexed_reference_fasta
      output_basename: output_basename
      tool_name: 
        valueFrom: ${return "vardict";}
      strip_info: strip_info
    out: [normalized_vcf]

  prep_mnp_variants:
    run: ../tools/prep_mnp_variants.cwl
    in:
      strelka2_vcf: normalize_strelka2_vcf/normalized_vcf
      other_vcfs: [normalize_mutect2_vcf/normalized_vcf, normalize_lancet_vcf/normalized_vcf, normalize_vardict_vcf/normalized_vcf]
      output_basename: output_basename
    out: [output_vcfs]

  bcbio_variant_recall_ensemble:
    run: ../tools/bcbio_variant_recall_ensemble.cwl
    in:
      input_vcfs: prep_mnp_variants/output_vcfs
      reference: indexed_reference_fasta
      min_overlap: min_overlap
      tool_name_csv: 
        valueFrom: ${return "strelka2,mutect2,lancet,vardict"}
      output_basename: output_basename
    out: [consensus_vcf]
