cwlVersion: v1.0
class: Workflow
id: kfdrc_somatic_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
inputs:
  reference_fasta: {type: File, doc: "Genome reference fasta file"}
  reference_dict: File
  wgs_calling_interval_list: File
  af_only_gnomad_vcf: {type: File, doc: "Broad GATK gnomad reference file"}
  exac_common_vcf: {type: File, doc: "Broad GATK exac reference file"}
  hg38_strelka_bed: {type: File, doc: "bgzipped chromosome bed file"}
  input_tumor_aligned:
    type: File
    secondaryFiles: |
      ${
        var dpath = self.location.replace(self.basename, "")
        if(self.nameext == '.bam'){
          return {"location": dpath+self.nameroot+".bai", "class": "File"}
        }
        else{
          return {"location": dpath+self.basename+".crai", "class": "File"}
        }
      }
    doc: "tumor BAM or CRAM"

  input_tumor_name: string
  input_normal_aligned:
    type: File
    secondaryFiles: |
      ${
        var dpath = self.location.replace(self.basename, "")
        if(self.nameext == '.bam'){
          return {"location": dpath+self.nameroot+".bai", "class": "File"}
        }
        else{
          return {"location": dpath+self.basename+".crai", "class": "File"}
        }
      }
    doc: "normal BAM or CRAM"

  input_normal_name: string
  threads: {type: int, doc: "For ControlFreeC.  Recommend 16 max, as I/O gets saturated after that losing any advantage."}
  # exome_flag: {type: ['null', string], doc: "insert 'Y' if exome mode"}
  # capture_regions: {type: ['null', File], doc: "If not WGS, provide this bed file"}
  select_vars_mode: {type: string, doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression"}
  vep_cache: {type: File, doc: "tar gzipped cache from ensembl/local converted cache" }
  chr_len: {type: File, doc: "file with chromosome lengths"}
  output_basename: string
  ploidy: {type: 'int[]', doc: "Array of ploidy possibilities for ControlFreeC to try"}
  mate_orientation_sample: {type: ['null', {type: enum, name: mate_orientation_sample, symbols: ["0", "FR", "RF", "FF"]}], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  mate_orientation_control: {type: ['null', {type: enum, name: mate_orientation_control, symbols: ["0", "FR", "RF", "FF"]}], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  b_allele: {type: ['null', File], doc: "germline calls, needed for BAF.  VarDict input recommended.  Tool will prefilter for germline and pass if expression given"}
  chr_len: {type: File, doc: "TSV with chromsome names and lengths. Limit to chromosome you actually want analyzed"}
  coeff_var: {type: float, default: 0.05, doc: "Coefficient of variantion to set window size.  Default 0.05 recommended"}
  contamination_adjustment: {type: ['null', boolean], doc: "TRUE or FALSE to have ControlFreec estimate normal contam"}
  include_expression: {type: ['null', string], doc: "Filter expression if vcf has mixed somatic/germline calls, use as-needed"}
  exclude_expression: {type: ['null', string], doc: "Filter expression if vcf has mixed somatic/germline calls, use as-needed"}
  sex: {type: ['null', {type: enum, name: sex, symbols: ["XX", "XY"] }], doc: "If known, XX for female, XY for male"}

outputs:
  strelka2_vep_vcf: {type: File, outputSource: vep_annot_strelka2/output_vcf}
  strelka2_vep_tbi: {type: File, outputSource: vep_annot_strelka2/output_tbi}
  strelka2_prepass_vcf: {type: File, outputSource: rename_strelka_samples/reheadered_vcf}
  strelka2_vep_maf: {type: File, outputSource: vep_annot_strelka2/output_maf}
  mutect2_vep_vcf: {type: File, outputSource: vep_annot_mutect2/output_vcf}
  mutect2_vep_tbi: {type: File, outputSource: vep_annot_mutect2/output_tbi}
  mutect2_prepass_vcf: {type: File, outputSource: filter_mutect2_vcf/filtered_vcf}
  mutect2_vep_maf: {type: File, outputSource: vep_annot_mutect2/output_maf}
  manta_vep_vcf: {type: File, outputSource: vep_annot_manta/output_vcf}
  manta_vep_tbi: {type: File, outputSource: vep_annot_manta/output_tbi}
  manta_prepass_vcf: {type: File, outputSource: rename_manta_samples/reheadered_vcf}
  manta_vep_maf: {type: File, outputSource: vep_annot_manta/output_maf}
  ctrlfreec_pval: {type: File, outputSource: rename_outputs/ctrlfreec_pval}
  ctrlfreec_config: {type: File, outputSource: rename_outputs/ctrlfreec_config}
  ctrlfreec_pngs: {type: 'File[]', outputSource: rename_outputs/ctrlfreec_pngs}
  ctrlfreec_bam_ratio: {type: File, outputSource: rename_outputs/ctrlfreec_bam_ratio}
  ctrlfreec_bam_seg: {type: File, outputSource: convert_ratio_to_seg/ctrlfreec_ratio2seg}
  ctrlfreec_baf: {type: File, outputSource: rename_outputs/ctrlfreec_baf}
  ctrlfreec_info: {type: File, outputSource: rename_outputs/ctrlfreec_info}

steps:

  index_references:
    run: ../tools/index_references.cwl
    in:
      input_fasta_file: reference_fasta
      af_only_gnomad_vcf: af_only_gnomad_vcf
      exac_common_vcf: exac_common_vcf
      hg38_strelka_bed: hg38_strelka_bed
    out:
      [indexed_reference_fasta, indexed_af_only_gnomad_vcf, indexed_exac_common_vcf, indexed_hg38_strelka_bed, reference_fai]


  bcftools_filter_vcf:
    run: ../tools/bcftools_filter_vcf.cwl
    in:
      input_vcf: b_allele
      include_expression: include_expression
      exclude_expression: exclude_expression
      output_basename: output_basename
    out:
      [filtered_vcf]

  controlfreec_tumor_mini_pileup:
    run: ../tools/control_freec_mini_pileup.cwl
    in:
      input_reads: samtools_tumor_cram2bam/bam_file
      threads:
        valueFrom: ${return 16}
      reference: index_references/indexed_reference_fasta
      snp_vcf: b_allele
    out:
      [pileup]

  controlfreec_normal_mini_pileup:
    run: ../tools/control_freec_mini_pileup.cwl
    in:
      input_reads: samtools_normal_cram2bam/bam_file
      threads:
        valueFrom: ${return 16}
      reference: index_references/indexed_reference_fasta
      snp_vcf: b_allele
    out:
      [pileup]

  samtools_tumor_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_tumor_aligned
      threads:
        valueFrom: ${return 36}
      reference: index_references/indexed_reference_fasta
    out: [bam_file]

  samtools_normal_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_normal_aligned
      threads:
        valueFrom: ${return 36}
      reference: index_references/indexed_reference_fasta
    out: [bam_file]

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
      # capture_regions: capture_regions
      max_threads: threads
      reference: index_references/indexed_reference_fasta
      snp_file: bcftools_filter_vcf/filtered_vcf
      coeff_var: coeff_var
      sex: sex
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
      reference_fai: index_references/reference_fai
      ctrlfreec_ratio: control_free_c/ratio
      sample_name: input_tumor_name
      output_basename: output_basename
    out: [ctrlfreec_ratio2seg]

  gatk_intervallisttools:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: wgs_calling_interval_list
      reference_dict: reference_dict
      # exome_flag: exome_flag
      scatter_ct:
        valueFrom: ${return 50}
      bands:
        valueFrom: ${return 80000000}
    out: [output]

  strelka2:
    run: ../tools/strelka2.cwl
    in:
      input_tumor_aligned: input_tumor_aligned
      input_normal_aligned: input_normal_aligned
      reference: index_references/indexed_reference_fasta
      hg38_strelka_bed: hg38_strelka_bed
      # exome_flag: exome_flag
    out: [output_snv, output_indel]

  manta:
    run: ../tools/manta.cwl
    in:
      input_tumor_cram: input_tumor_aligned
      input_normal_cram: input_normal_aligned
      output_basename: output_basename
      reference: index_references/indexed_reference_fasta
      hg38_strelka_bed: hg38_strelka_bed
    out: [output_sv]
  
  mutect2:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge;ebs-gp2;500
    run: ../tools/gatk_Mutect2.cwl
    in:
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      reference: index_references/indexed_reference_fasta
      reference_dict: reference_dict
      interval_list: gatk_intervallisttools/output
      af_only_gnomad_vcf: index_references/indexed_af_only_gnomad_vcf
      # exome_flag: exome_flag
    scatter: [interval_list]
    out: [mutect2_vcf, f1r2_counts, mutect_stats]
  
  mutect2_filter_support:
    run: ../workflow/kfdrc_mutect2_filter_support_subwf.cwl
    in:
      indexed_reference_fasta: index_references/indexed_reference_fasta
      reference_dict: reference_dict
      wgs_calling_interval_list: gatk_intervallisttools/output
      input_tumor_aligned: input_tumor_aligned
      input_normal_aligned: input_normal_aligned
      exac_common_vcf: index_references/indexed_exac_common_vcf
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

  rename_manta_samples:
    run: ../tools/bcftools_reheader_vcf.cwl
    in:
      input_vcf: manta/output_sv
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
      reference: index_references/indexed_reference_fasta
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

  gatk_selectvariants_manta:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    run: ../tools/gatk_selectvariants.cwl
    label: GATK Select Manta PASS
    in:
      input_vcf: rename_manta_samples/reheadered_vcf
      output_basename: output_basename
      tool_name:
        valueFrom: ${return "manta"}
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
      reference: index_references/indexed_reference_fasta
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
      reference: index_references/indexed_reference_fasta
      cache: vep_cache
    out: [output_vcf, output_tbi, output_maf, warn_txt]

  vep_annot_manta:
    run: ../tools/vep_vcf2maf.cwl
    in:
      input_vcf: gatk_selectvariants_manta/pass_vcf
      output_basename: output_basename
      tumor_id: input_tumor_name
      normal_id: input_normal_name
      tool_name:
        valueFrom: ${return "manta_somatic"}
      reference: index_references/indexed_reference_fasta
      cache: vep_cache
    out: [output_vcf, output_tbi, output_maf, warn_txt]
  
$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4