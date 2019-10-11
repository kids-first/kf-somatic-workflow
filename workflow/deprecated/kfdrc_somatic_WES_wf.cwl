cwlVersion: v1.0
class: Workflow
id: kfdrc_whole_exome_somatic_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict]}
  reference_dict: File
  threads: {type: int, doc: "For ControlFreeC.  Recommend 16 max, as I/O gets saturated after that losing any advantage."}
  chr_len: {type: File, doc: "Text file with lengths of each chromosome in fasta file"}
  ref_chrs: {type: File, doc: "Tar gzipped file, one fasta file per chromosome. For now, hardcoded to create folder GRCh38_everyChrs when unpacked"}
  calling_interval_list: {type: File, doc: "Bed preferred for best compatibility Mutect2, Strelka2, etc"}
  af_only_gnomad_vcf: {type: File, secondaryFiles: ['.tbi']}
  exac_common_vcf: {type: File, secondaryFiles: ['.tbi']}
  capture_regions: {type: File, doc: "Bed file for CNV calls (exact intervals)"}
  b_allele_vcf: {type: File, doc: "vcf containing SNV b-alleles sites (only sites with PASS will be used)"}
  hg38_strelka_bed: File
  manifest: {type: File, doc: "Nextera Manifest file for Canvas"}
  sample_name: string
  canvas_reference_file: {type: File, doc: "Canvas-ready kmer file"} 
  genomeSize_file: {type: File, doc: "GenomeSize.xml file"}
  genome_fasta: {type: File, doc: "Genome.fa", secondaryFiles: [.fai]}
  input_tumor_aligned: { type: File, secondaryFiles: [.crai] }
  input_tumor_name: string
  input_normal_aligned: { type: File, secondaryFiles: [.crai] }
  input_normal_name: string
  exome_flag: string
  select_vars_mode: {type: string, doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression"}
  vep_cache: {type: File, label: tar gzipped cache from ensembl/local converted cache}
  output_basename: string

outputs:
  strelka2_vep_vcf: {type: File, outputSource: vep_annot_strelka2/output_vcf}
  strelka2_vep_tbi: {type: File, outputSource: vep_annot_strelka2/output_tbi}
  strelka2_prepass_vcf: {type: File, outputSource: rename_strelka_samples/reheadered_vcf}
  strelka2_vep_maf: {type: File, outputSource: vep_annot_strelka2/output_maf}
  mutect2_vep_vcf: {type: File, outputSource: vep_annot_mutect2/output_vcf}
  mutect2_vep_tbi: {type: File, outputSource: vep_annot_mutect2/output_tbi}
  mutect2_prepass_vcf: {type: File, outputSource: filter_mutect2_vcf/filtered_vcf}
  mutect2_vep_maf: {type: File, outputSource: vep_annot_mutect2/output_maf}
  ctrlfreec_cnv: {type: File, outputSource: control_free_c/output_cnv}
  ctrlfreec_cnv_bam_ratio: { type: File, outputSource: control_free_c/output_txt }
  ctrlfreec_cnv_pval: { type: File, outputSource: control_free_c_r/output_pval }
  ctrlfreec_cnv_png: { type: File, outputSource: control_free_c_viz/output_png }
  canvas_cnv_vcf: {type: File, outputSource: canvas/output_vcf}
  canvas_coverage_txt: {type: File, outputSource: canvas/output_txt}
  canvas_folder: {type: File, outputSource: canvas/output_folder}

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

  gen_config:
    run: ../tools/gen_controlfreec_configfile.cwl
    in:
      tumor_bam: samtools_tumor_cram2bam/bam_file
      normal_bam: samtools_normal_cram2bam/bam_file 
      reference: indexed_reference_fasta
      capture_regions: capture_regions
      exome_flag: exome_flag
      chr_len: chr_len
      threads: threads
    out: [config_file]

  control_free_c:
    run: ../tools/control_freec.cwl
    in:
      tumor_bam: samtools_tumor_cram2bam/bam_file
      normal_bam: samtools_normal_cram2bam/bam_file
      ref_chrs: ref_chrs
      chr_len: chr_len
      threads: threads
      config_file: gen_config/config_file
      capture_regions: capture_regions
      output_basename: output_basename
    out: [output_txt, output_cnv]
  
  control_free_c_r:
    run: ../tools/control_freec_R.cwl
    in:
      cnv_bam_ratio: control_free_c/output_txt
      cnv_result: control_free_c/output_cnv
    out: [output]

  control_free_c_viz:
    run: ../tools/control_freec_visualize.cwl
    in:
      output_basename: output_basename
      cnv_bam_ratio: control_free_c/output_txt
    out: [output]

  canvas:
    run: ../tools/canvas-paired-wes.cwl
    in:  
      tumor_bam: samtools_tumor_cram2bam/bam_file
      control_bam: samtools_normal_cram2bam/bam_file
      manifest: manifest
      b_allele_vcf: b_allele_vcf
      sample_name: sample_name
      output_basename: output_basename
      reference: canvas_reference_file
      genomeSize_file: genomeSize_file
      genome_fasta: genome_fasta
      filter_bed: capture_regions
    out: [output_vcf, output_txt, output_folder]

  gatk_intervallisttools:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: calling_interval_list
      reference_dict: reference_dict
      exome_flag: exome_flag
      scatter_ct:
        valueFrom: ${return 50}
      bands:
        valueFrom: ${return 80000000}
    out: [output]

  strelka2:
    run: ../tools/strelka2.cwl
    in:
      input_tumor_aligned: samtools_tumor_cram2bam/bam_file
      input_normal_aligned: samtools_normal_cram2bam/bam_file
      reference: indexed_reference_fasta
      hg38_strelka_bed: hg38_strelka_bed
      exome_flag: exome_flag
    out: [output]
  
  mutect2:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge;ebs-gp2;500
    run: ../tools/gatk_Mutect2.cwl
    in:
      input_tumor_aligned: samtools_tumor_cram2bam/bam_file
      input_tumor_name: input_tumor_name
      input_normal_aligned: samtools_normal_cram2bam/bam_file
      input_normal_name: input_normal_name
      reference: indexed_reference_fasta
      interval_list: gatk_intervallisttools/output
      af_only_gnomad_vcf: af_only_gnomad_vcf
      exome_flag: exome_flag
    scatter: [interval_list]
    out: [mutect2_vcf, f1r2_counts, mutect_stats]
  
  mutect2_filter_support:
    run: ../workflow/kfdrc_mutect2_filter_support_subwf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      reference_dict: reference_dict
      wgs_calling_interval_list: gatk_intervallisttools/output
      input_tumor_aligned: samtools_tumor_cram2bam/bam_file
      input_normal_aligned: samtools_normal_cram2bam/bam_file
      exac_common_vcf: exac_common_vcf
      output_basename: output_basename
      f1r2_counts: mutect2/f1r2_counts
    out: [contamination_table, segmentation_table, f1r2_bias]
  
  merge_strelka2_vcf:
    run: ../tools/gatk_mergevcfs.cwl
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
    run: ../tools/gatk_mergevcfs.cwl
    label: Merge & pass filter mutect2
    in:
      input_vcfs: mutect2/mutect2_vcf
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name:
        valueFrom: ${return "mutect2"}
    out: [merged_vcf]

  merge_mutect2_stats:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    run: ../tools/gatk_mergemutectstats.cwl
    label: Merge mutect2 stats
    in:
      input_stats: mutect2/mutect_stats
      output_basename: output_basename
    out: [merged_stats]
  
  filter_mutect2_vcf:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    run: ../tools/gatk_filtermutectcalls.cwl
    in:
      mutect_vcf: merge_mutect2_vcf/merged_vcf
      mutect_stats: merge_mutect2_stats/merged_stats
      reference: indexed_reference_fasta
      output_basename: output_basename
      contamination_table: mutect2_filter_support/contamination_table
      segmentation_table: mutect2_filter_support/segmentation_table
      ob_priors: mutect2_filter_support/f1r2_bias
    out: [stats_table, filtered_vcf]

  gatk_selectvariants_mutect2:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    run: ../tools/gatk_selectvariants.cwl
    label: GATK Select Mutect2 PASS
    in:
      input_vcf: filter_mutect2_vcf/filtered_vcf
      output_basename: output_basename
      tool_name:
        valueFrom: ${return "mutect2"}
      mode: select_vars_mode
    out: [pass_vcf]

  gatk_selectvariants_strelka2:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    run: ../tools/gatk_selectvariants.cwl
    label: GATK Select Strelka2 PASS
    in:
      input_vcf: rename_strelka_samples/reheadered_vcf
      output_basename: output_basename
      tool_name:
        valueFrom: ${return "strelka2"}
      mode: select_vars_mode
    out: [pass_vcf]
    
  vep_annot_strelka2:
    run: ../tools/vep_vcf2maf.cwl
    in:
      input_vcf: gatk_selectvariants_strelka2/pass_vcf
      output_basename: output_basename
      tumor_id: input_tumor_name
      normal_id: input_normal_name
      tool_name:
        valueFrom: ${return "strelka2_somatic"}
      reference: indexed_reference_fasta
      cache: vep_cache
    out: [output_vcf, output_tbi, output_maf, warn_txt]

  vep_annot_mutect2:
    run: ../tools/vep_vcf2maf.cwl
    in:
      input_vcf: gatk_selectvariants_mutect2/pass_vcf
      output_basename: output_basename
      tumor_id: input_tumor_name
      normal_id: input_normal_name
      tool_name:
        valueFrom: ${return "mutect2_somatic"}
      reference: indexed_reference_fasta
      cache: vep_cache
    out: [output_vcf, output_tbi, output_maf, warn_txt]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4