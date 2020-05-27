cwlVersion: v1.0
class: Workflow
id: kfdrc_prepare_reference
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  input_fasta: File
  input_fai: { type: 'File?' }
  input_dict: { type: 'File?' }
  input_alt: { type: 'File?' }
  input_amb: { type: 'File?' }
  input_ann: { type: 'File?' }
  input_bwt: { type: 'File?' }
  input_pac: { type: 'File?' }
  input_sa: { type: 'File?' }
  generate_bwa_indexes: { type: 'boolean?' }

outputs:
  indexed_fasta: { type: File, doc: "Reference fasta with all available indexes as secondaryFiles", outputSource: bundle_secondaries/output } 
  reference_dict: { type: File, doc: "Standalone reference dict", outputSource: picard_create_sequence_dictionary/dict }

steps:
  samtools_faidx:
    run: ../tools/samtools_faidx.cwl
    in:
      input_fasta: input_fasta
      input_index: input_fai
    out: [fai]
  picard_create_sequence_dictionary:
    run: ../tools/picard_createsequencedictionary.cwl
    in:
      input_fasta: input_fasta
      input_dict: input_dict
    out: [dict]
  bwa_index:
    run: ../tools/bwa_index.cwl
    in:
      generate_bwa_indexes: generate_bwa_indexes
      input_fasta: input_fasta
      input_alt: input_alt
      input_amb: input_amb
      input_ann: input_ann
      input_bwt: input_bwt
      input_pac: input_pac
      input_sa: input_sa
    out: [alt,amb,ann,bwt,pac,sa]
  bundle_secondaries:
    run: ../tools/bundle_secondaryfiles.cwl
    in:
      primary_file: input_fasta
      secondary_files:
        source: [samtools_faidx/fai,picard_create_sequence_dictionary/dict,bwa_index/alt,bwa_index/amb,bwa_index/ann,bwa_index/bwt,bwa_index/pac,bwa_index/sa]
        linkMerge: merge_flattened
    out: [output]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: sbg:maxNumberOfParallelInstances
    value: 2
