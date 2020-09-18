cwlVersion: v1.0
class: Workflow
id: kfdrc_controlfreec_wf

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  input_tumor: {type: File}
  input_normal: {type: File}
  input_tumor_name: {type: string, doc: "Sample name to put into the converted seg file"}
  threads: {type: int, doc: "Number of threads to run controlfreec.  Going above 16 is not recommended, there is no apparent added value"}
  output_basename: string
  ploidy: {type: 'int[]', doc: "Array of ploidy possibilities for ControlFreeC to try"}
  mate_orientation_sample: {type: ['null', {type: enum, name: mate_orientation_sample, symbols: ["0", "FR", "RF", "FF"]}], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  mate_orientation_control: {type: ['null', {type: enum, name: mate_orientation_control, symbols: ["0", "FR", "RF", "FF"]}], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  capture_regions: {type: ['null', File], doc: "If not WGS, provide "}
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai]}
  reference_fai: {type: File, doc: "fasta index file for seg file conversion"}
  b_allele: {type: ['null', File], doc: "germline calls, needed for BAF.  GATK HC VQSR input recommended.  Tool will prefilter for germline and pass if expression given"}
  chr_len: {type: File, doc: "TSV with chromsome names and lengths. Limit to chromosome you actually want analyzed"}
  coeff_var: {type: float, default: 0.05, doc: "Coefficient of variantion to set window size.  Default 0.05 recommended"}
  contamination_adjustment: {type: ['null', boolean], doc: "TRUE or FALSE to have ControlFreec estimate normal contam"}
  cfree_sex: {type: ['null', {type: enum, name: sex, symbols: ["XX", "XY"] }], doc: "If known, XX for female, XY for male"}

outputs:
  ctrlfreec_pval: {type: File, outputSource: rename_outputs/ctrlfreec_pval}
  ctrlfreec_config: {type: File, outputSource: rename_outputs/ctrlfreec_config}
  ctrlfreec_pngs: {type: 'File[]', outputSource: rename_outputs/ctrlfreec_pngs}
  ctrlfreec_bam_ratio: {type: File, outputSource: rename_outputs/ctrlfreec_bam_ratio}
  ctrlfreec_bam_seg: {type: File, outputSource: convert_ratio_to_seg/ctrlfreec_ratio2seg}
  ctrlfreec_baf: {type: File, outputSource: rename_outputs/ctrlfreec_baf}
  ctrlfreec_info: {type: File, outputSource: rename_outputs/ctrlfreec_info}

steps:
  gatk_filter_germline:
    run: ../tools/gatk_filter_germline_variant.cwl
    in:
      input_vcf: b_allele
      reference_fasta: indexed_reference_fasta
      output_basename: output_basename
    out:
      [filtered_vcf, filtered_pass_vcf]

  controlfreec_tumor_mini_pileup:
    run: ../tools/control_freec_mini_pileup.cwl
    in:
      input_reads: samtools_tumor_cram2bam/bam_file
      threads:
        valueFrom: ${return 16}
      reference: indexed_reference_fasta
      snp_vcf: b_allele
    out:
      [pileup]

  controlfreec_normal_mini_pileup:
    run: ../tools/control_freec_mini_pileup.cwl
    in:
      input_reads: samtools_normal_cram2bam/bam_file
      threads:
        valueFrom: ${return 16}
      reference: indexed_reference_fasta
      snp_vcf: b_allele
    out:
      [pileup]

  samtools_tumor_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_tumor
      threads:
        valueFrom: ${return 16}
      reference: indexed_reference_fasta
    out:
      [bam_file]

  samtools_normal_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_normal
      threads:
        valueFrom: ${return 16}
      reference: indexed_reference_fasta
    out:
      [bam_file]

  control_free_c: 
    run: ../tools/control-freec-11-6-sbg.cwl
    in: 
      mate_file_sample: samtools_tumor_cram2bam/bam_file
      mate_orientation_sample: mate_orientation_sample
      mini_pileup_sample: controlfreec_tumor_mini_pileup/pileup
      mate_file_control: samtools_normal_cram2bam/bam_file
      mate_orientation_control: mate_orientation_control
      mini_pileup_control: controlfreec_normal_mini_pileup/pileup
      chr_len: chr_len
      ploidy: ploidy
      capture_regions: capture_regions
      max_threads: threads
      reference: indexed_reference_fasta
      snp_file: gatk_filter_germline/filtered_pass_vcf
      coeff_var: coeff_var
      sex: cfree_sex
      contamination_adjustment: contamination_adjustment
    out: [cnvs_pvalue, config_script, pngs, ratio, sample_BAF, info_txt]

  rename_outputs:
    run: ../tools/ubuntu_rename_outputs.cwl
    in:
      input_files: [control_free_c/cnvs_pvalue, control_free_c/config_script, control_free_c/ratio, control_free_c/sample_BAF, control_free_c/info_txt]
      input_pngs: control_free_c/pngs
      output_basename: output_basename
    out: [ctrlfreec_pval, ctrlfreec_config, ctrlfreec_pngs, ctrlfreec_bam_ratio, ctrlfreec_baf, ctrlfreec_info]
  
  convert_ratio_to_seg:
    run: ../tools/ubuntu_ratio2seg.cwl
    in:
      reference_fai: reference_fai
      ctrlfreec_ratio: control_free_c/ratio
      sample_name: input_tumor_name
      output_basename: output_basename
    out: [ctrlfreec_ratio2seg]
$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 2