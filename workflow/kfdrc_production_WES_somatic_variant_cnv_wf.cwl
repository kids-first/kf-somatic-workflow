cwlVersion: v1.0
class: Workflow
id: kfdrc_production_somatic_wes_variant_cnv_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict]}
  reference_fai: File
  reference_dict: File
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
  exome_flag: {type: string?, default: "Y", doc: "Whether to run in exome mode for callers. Should be Y or leave blank as default is Y. Only make N if you are certain"}
  cfree_threads: {type: ['null', int], doc: "For ControlFreeC.  Recommend 16 max, as I/O gets saturated after that losing any advantage.", default: 16}
  vardict_min_vaf: {type: ['null', float], doc: "Min variant allele frequency for vardict to consider.  Recommend 0.05", default: 0.05}
  select_vars_mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression", default: "gatk"}
  vardict_cpus: {type: ['null', int], default: 9}
  vardict_ram: {type: ['null', int], default: 18, doc: "In GB"}
  lancet_ram: {type: ['null', int], default: 12, doc: "Adjust in rare circumstances in which 12 GB is not enough"}
  lancet_window: {type: ['null', int], doc: "window size for lancet.  default is 600, recommend 500 for WGS, 600 for exome+", default: 600}
  lancet_padding: {type: ['null', int], doc: "Recommend 0 if interval file padded already, half window size if not", default: 0}
  vardict_padding: {type: ['null', int], doc: "Padding to add to input intervals, recommend 0 if intervals already padded, 150 if not", default: 0}
  vep_cache: {type: File, doc: "tar gzipped cache from ensembl/local converted cache"}
  vep_ref_build: {type: ['null', string], doc: "Genome ref build used, should line up with cache.", default: "GRCh38" }
  padded_capture_regions: {type: ['null', File], doc: "Recommend 100bp pad, for somatic variant"}
  unpadded_capture_regions: {type: ['null', File], doc: "Capture regions with NO padding for cnv calling"}
  cfree_chr_len: {type: File, doc: "file with chromosome lengths"}
  cfree_ploidy: {type: 'int[]', doc: "Array of ploidy possibilities for ControlFreeC to try"}
  cfree_mate_orientation_sample: {type: ['null', {type: enum, name: mate_orientation_sample, symbols: ["0", "FR", "RF", "FF"]}], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  cfree_mate_orientation_control: {type: ['null', {type: enum, name: mate_orientation_control, symbols: ["0", "FR", "RF", "FF"]}], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  b_allele: {type: ['null', File], secondaryFiles: ['.tbi'],  doc: "germline calls, needed for BAF.  GATK HC VQSR input recommended.  Tool will prefilter for germline and pass if expression given"}
  i_flag: {type: ['null', string], doc: "Flag to intersect germline calls on padded regions.  Use N if you want to skip this."}
  mutect2_af_only_gnomad_vcf: {type: File, secondaryFiles: ['.tbi']}
  mutect2_exac_common_vcf: {type: File, secondaryFiles: ['.tbi']}
  hg38_strelka_bed: {type: File, secondaryFiles: ['.tbi'], doc: "Bgzipped interval bed file. Recommned padding 100bp"}
  cnvkit_annotation_file: {type: File, doc: "refFlat.txt file"}
  cfree_coeff_var: {type: ['null', float], default: 0.05, doc: "Coefficient of variation to set window size.  Default 0.05 recommended"}
  cfree_contamination_adjustment: {type: ['null', boolean], doc: "TRUE or FALSE to have ControlFreec estimate normal contam"}
  # combined_include_expression: {type: ['null', string], doc: "Filter expression if vcf has non-PASS combined calls, use as-needed"}
  # combined_exclude_expression: {type: ['null', string], doc: "Filter expression if vcf has non-PASS combined calls, use as-needed"}  
  # min_theta2_frac: {type: ['null', float], doc: "Minimum fraction of genome with copy umber alterations.  Default is 0.05, recommend 0.01", default: 0.01}
  cfree_sex: {type: ['null', {type: enum, name: sex, symbols: ["XX", "XY"] }], doc: "If known, XX for female, XY for male"}
  cnvkit_sex: {type: ['null', string ], doc: "If known, choices are m,y,male,Male,f,x,female,Female"}
  output_basename: string

outputs:
  ctrlfreec_pval: {type: File, outputSource: run_controlfreec/ctrlfreec_pval}
  ctrlfreec_config: {type: File, outputSource: run_controlfreec/ctrlfreec_config}
  ctrlfreec_pngs: {type: 'File[]', outputSource: run_controlfreec/ctrlfreec_pngs}
  ctrlfreec_bam_ratio: {type: File, outputSource: run_controlfreec/ctrlfreec_bam_ratio}
  ctrlfreec_bam_seg: {type: File, outputSource: run_controlfreec/ctrlfreec_bam_seg}
  ctrlfreec_baf: {type: File, outputSource: run_controlfreec/ctrlfreec_baf}
  ctrlfreec_info: {type: File, outputSource: run_controlfreec/ctrlfreec_info}
  cnvkit_cnr: {type: File, outputSource: run_cnvkit/cnvkit_cnr}
  cnvkit_cnn_output: {type: ['null', File], outputSource: run_cnvkit/cnvkit_cnn_output}
  cnvkit_calls: {type: File, outputSource: run_cnvkit/cnvkit_calls}
  cnvkit_metrics: {type: File, outputSource: run_cnvkit/cnvkit_metrics}
  cnvkit_gainloss: {type: File, outputSource: run_cnvkit/cnvkit_gainloss}
  cnvkit_seg: {type: File, outputSource: run_cnvkit/cnvkit_seg}
#   theta2_calls: {type: File, outputSource: cnvkit_import_theta2/theta2_adjusted_cns}
#   theta2_seg: {type: File, outputSource: cnvkit_import_theta2/theta2_adjusted_seg}
#   theta2_subclonal_results: {type: 'File[]', outputSource: [run_theta2/n3_graph, run_theta2/n2_results, run_theta2/best_results]}
#   theta2_subclonal_cns: {type: 'File[]', outputSource: cnvkit_import_theta2/theta2_subclone_cns}
#   theta2_subclone_seg: {type: 'File[]', outputSource: cnvkit_import_theta2/theta2_subclone_seg}
  strelka2_vep_vcf: {type: File, outputSource: run_strelka2/strelka2_vep_vcf}
  strelka2_vep_tbi: {type: File, outputSource: run_strelka2/strelka2_vep_tbi}
  strelka2_prepass_vcf: {type: File, outputSource: run_strelka2/strelka2_prepass_vcf}
  strelka2_vep_maf: {type: File, outputSource: run_strelka2/strelka2_vep_maf}
  mutect2_vep_vcf: {type: File, outputSource: run_mutect2/mutect2_vep_vcf}
  mutect2_vep_tbi: {type: File, outputSource: run_mutect2/mutect2_vep_tbi}
  mutect2_prepass_vcf: {type: File, outputSource: run_mutect2/mutect2_filtered_vcf}
  mutect2_vep_maf: {type: File, outputSource: run_mutect2/mutect2_vep_maf}
  vardict_vep_somatic_only_vcf: {type: File, outputSource: run_vardict/vardict_vep_somatic_only_vcf}
  vardict_vep_somatic_only_tbi: {type: File, outputSource: run_vardict/vardict_vep_somatic_only_tbi}
  vardict_vep_somatic_only_maf: {type: File, outputSource: run_vardict/vardict_vep_somatic_only_maf}
  vardict_prepass_vcf: {type: File, outputSource: run_vardict/vardict_prepass_vcf}
  lancet_vep_vcf: {type: File, outputSource: run_lancet/lancet_vep_vcf}
  lancet_vep_tbi: {type: File, outputSource: run_lancet/lancet_vep_tbi}
  lancet_vep_maf: {type: File, outputSource: run_lancet/lancet_vep_maf}
  lancet_prepass_vcf: {type: File, outputSource: run_lancet/lancet_prepass_vcf}

steps:
  bedtools_intersect_germline:
    run: ../tools/bedtools_intersect.cwl
    in:
      input_vcf: b_allele
      output_basename: output_basename
      input_bed_file: unpadded_capture_regions
      flag: i_flag
    out:
      [intersected_vcf]

  gatk_intervallisttools:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: padded_capture_regions
      reference_dict: reference_dict
      exome_flag: exome_flag
      scatter_ct:
        valueFrom: ${return 50}
      bands:
        valueFrom: ${return 80000000}
    out: [output]
  
  gatk_filter_germline:
    run: ../tools/gatk_filter_germline_variant.cwl
    in:
      input_vcf: bedtools_intersect_germline/intersected_vcf
      reference_fasta: indexed_reference_fasta
      output_basename: output_basename
    out:
      [filtered_vcf, filtered_pass_vcf]

  samtools_cram2bam_plus_calmd_tumor:
    run: ../tools/samtools_cram2bam_plus_calmd.cwl
    in:
      input_reads: input_tumor_aligned
      threads:
        valueFrom: ${return 16;}
      reference: indexed_reference_fasta
    out: [bam_file]

  samtools_cram2bam_plus_calmd_normal:
    run: ../tools/samtools_cram2bam_plus_calmd.cwl
    in:
      input_reads: input_normal_aligned
      threads:
        valueFrom: ${return 16;}
      reference: indexed_reference_fasta
    out: [bam_file]

  run_vardict:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge
    run: ../sub_workflows/kfdrc_vardict_sub_wf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_tumor_aligned: samtools_cram2bam_plus_calmd_tumor/bam_file
      input_tumor_name: input_tumor_name
      input_normal_aligned: samtools_cram2bam_plus_calmd_normal/bam_file
      input_normal_name: input_normal_name
      output_basename: output_basename
      reference_dict: reference_dict
      bed_invtl_split: gatk_intervallisttools/output
      padding: vardict_padding
      min_vaf: vardict_min_vaf
      select_vars_mode: select_vars_mode
      cpus: vardict_cpus
      ram: vardict_ram
      vep_cache: vep_cache
      vep_ref_build: vep_ref_build
    out:
      [vardict_vep_somatic_only_vcf, vardict_vep_somatic_only_tbi, vardict_vep_somatic_only_maf, vardict_prepass_vcf]

  run_mutect2:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge
    run: ../sub_workflows/kfdrc_mutect2_sub_wf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      reference_dict: reference_dict
      bed_invtl_split: gatk_intervallisttools/output
      af_only_gnomad_vcf: mutect2_af_only_gnomad_vcf
      exac_common_vcf: mutect2_exac_common_vcf
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      exome_flag: exome_flag
      vep_cache: vep_cache
      vep_ref_build: vep_ref_build
      output_basename: output_basename
      select_vars_mode: select_vars_mode
    out:
      [mutect2_filtered_stats, mutect2_filtered_vcf, mutect2_vep_vcf, mutect2_vep_tbi, mutect2_vep_maf]

  run_strelka2:
    run: ../sub_workflows/kfdrc_strelka2_sub_wf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      reference_dict: reference_dict
      hg38_strelka_bed: hg38_strelka_bed
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      exome_flag: exome_flag
      vep_cache: vep_cache
      vep_ref_build: vep_ref_build
      output_basename: output_basename
      select_vars_mode: select_vars_mode
    out:
      [strelka2_vep_vcf, strelka2_vep_tbi, strelka2_prepass_vcf, strelka2_vep_maf]

  run_lancet:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge
    run: ../sub_workflows/kfdrc_lancet_sub_wf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_tumor_aligned: samtools_cram2bam_plus_calmd_tumor/bam_file
      input_tumor_name: input_tumor_name
      input_normal_aligned: samtools_cram2bam_plus_calmd_normal/bam_file
      input_normal_name: input_normal_name
      output_basename: output_basename
      select_vars_mode: select_vars_mode
      reference_dict: reference_dict
      bed_invtl_split: gatk_intervallisttools/output
      ram: lancet_ram
      window: lancet_window
      padding: lancet_padding
      vep_cache: vep_cache
      vep_ref_build: vep_ref_build
    out:
      [lancet_vep_vcf, lancet_vep_tbi, lancet_vep_maf, lancet_prepass_vcf]

  run_controlfreec:
    run: ../sub_workflows/kfdrc_controlfreec_sub_wf.cwl
    in:
      input_tumor_aligned: samtools_cram2bam_plus_calmd_tumor/bam_file
      input_tumor_name: input_tumor_name
      input_normal_aligned: samtools_cram2bam_plus_calmd_normal/bam_file
      threads: cfree_threads
      output_basename: output_basename
      ploidy: cfree_ploidy
      mate_orientation_sample: cfree_mate_orientation_sample
      mate_orientation_control: cfree_mate_orientation_control
      capture_regions: unpadded_capture_regions
      indexed_reference_fasta: indexed_reference_fasta
      reference_fai: reference_fai
      b_allele: gatk_filter_germline/filtered_pass_vcf
      chr_len: cfree_chr_len
      coeff_var: cfree_coeff_var
      contamination_adjustment: cfree_contamination_adjustment
      cfree_sex: cfree_sex
    out:
      [ctrlfreec_pval, ctrlfreec_config, ctrlfreec_pngs, ctrlfreec_bam_ratio, ctrlfreec_bam_seg, ctrlfreec_baf, ctrlfreec_info]

  run_cnvkit:
    run: ../sub_workflows/kfdrc_cnvkit_sub_wf.cwl
    in:
      input_tumor_aligned: samtools_cram2bam_plus_calmd_tumor/bam_file
      tumor_sample_name: input_tumor_name
      input_normal_aligned: samtools_cram2bam_plus_calmd_normal/bam_file
      reference: indexed_reference_fasta
      normal_sample_name: input_normal_name
      b_allele_vcf: gatk_filter_germline/filtered_pass_vcf
      capture_regions: unpadded_capture_regions
      annotation_file: cnvkit_annotation_file
      output_basename: output_basename
      sex: cnvkit_sex
    out:
      [cnvkit_cnr, cnvkit_cnn_output, cnvkit_calls, cnvkit_metrics, cnvkit_gainloss, cnvkit_seg]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4
