cwlVersion: v1.0
class: Workflow
id: kfdrc_production_somatic_wgs_variant_cnv_wf
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
  exome_flag: {type: string?, default: "N", doc: "Whether to run in exome mode for callers. Should be N or leave blank as default is N. Only make Y if you are certain"}
  wgs_calling_interval_list: {type: File, doc: "GATK intervals list-style, or bed file.  Recommend canocical chromosomes with N regions removed"}
  lancet_calling_interval_bed: {type: File, doc: "For WGS, highly recommended to use CDS bed, and supplement with region calls from strelka2 & mutect2.  Can still give calling list as bed if true WGS calling desired instead of exome+."}
  cfree_threads: {type: ['null', int], doc: "For ControlFreeC.  Recommend 16 max, as I/O gets saturated after that losing any advantage.", default: 16}
  vardict_min_vaf: {type: ['null', float], doc: "Min variant allele frequency for vardict to consider.  Recommend 0.05", default: 0.05}
  select_vars_mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression", default: "gatk"}
  vardict_cpus: {type: ['null', int], default: 9}
  vardict_ram: {type: ['null', int], default: 18, doc: "In GB"}
  lancet_ram: {type: ['null', int], default: 12, doc: "Adjust in rare circumstances in which 12 GB is not enough"}
  lancet_window: {type: ['null', int], doc: "window size for lancet.  default is 600, recommend 500 for WGS, 600 for exome+", default: 600}
  lancet_padding: {type: ['null', int], doc: "Recommend 0 if interval file padded already, half window size if not", default: 300}
  vardict_padding: {type: ['null', int], doc: "Padding to add to input intervals, recommened 0 if intervals already padded, 150 if not", default: 150}
  vep_cache: {type: File, doc: "tar gzipped cache from ensembl/local converted cache"}
  vep_ref_build: {type: ['null', string], doc: "Genome ref build used, should line up with cache.", default: "GRCh38" }
  cfree_chr_len: {type: File, doc: "file with chromosome lengths"}
  cfree_ploidy: {type: 'int[]', doc: "Array of ploidy possibilities for ControlFreeC to try"}
  cfree_mate_orientation_sample: {type: ['null', {type: enum, name: mate_orientation_sample, symbols: ["0", "FR", "RF", "FF"]}], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  cfree_mate_orientation_control: {type: ['null', {type: enum, name: mate_orientation_control, symbols: ["0", "FR", "RF", "FF"]}], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  b_allele: {type: ['null', File], secondaryFiles: ['.tbi'],  doc: "germline calls, needed for BAF.  GATK HC VQSR input recommended.  Tool will prefilter for germline and pass if expression given"}
  mutect2_af_only_gnomad_vcf: {type: File, secondaryFiles: ['.tbi']}
  mutect2_exac_common_vcf: {type: File, secondaryFiles: ['.tbi']}
  hg38_strelka_bed: {type: File, secondaryFiles: ['.tbi'], doc: "Bgzipped interval bed file. Recommend canonical chromosomes"}
  cnvkit_annotation_file: {type: File, doc: "refFlat.txt file"}
  cnvkit_wgs_mode: {type: ['null', string], doc: "for WGS mode, input Y. leave blank for hybrid mode", default: "Y"}
  cfree_coeff_var: {type: ['null', float], default: 0.05, doc: "Coefficient of variation to set window size.  Default 0.05 recommended"}
  cfree_contamination_adjustment: {type: ['null', boolean], doc: "TRUE or FALSE to have ControlFreec estimate normal contam"}
  combined_include_expression: {type: ['null', string], doc: "Filter expression if vcf has non-PASS combined calls, use as-needed, i.e. for VarDict: FILTER=\"PASS\" && (INFO/STATUS=\"Germline\" | INFO/STATUS=\"StrongSomatic\")"}
  combined_exclude_expression: {type: ['null', string], doc: "Filter expression if vcf has non-PASS combined calls, use as-needed"}  
  min_theta2_frac: {type: ['null', float], doc: "Minimum fraction of genome with copy umber alterations.  Default is 0.05, recommend 0.01", default: 0.01}
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
  theta2_calls: {type: File?, outputSource: run_theta2_purity/theta2_adjusted_cns}
  theta2_seg: {type: File?, outputSource: run_theta2_purity/theta2_adjusted_seg}
  theta2_subclonal_results: {type: ['null', 'File[]'], outputSource: run_theta2_purity/theta2_subclonal_results}
  theta2_subclonal_cns: {type: ['null', 'File[]'], outputSource: run_theta2_purity/theta2_subclonal_cns}
  theta2_subclone_seg: {type: ['null', 'File[]'], outputSource: run_theta2_purity/theta2_subclone_seg}
  strelka2_vep_vcf: {type: File, outputSource: run_strelka2/strelka2_vep_vcf}
  strelka2_vep_tbi: {type: File, outputSource: run_strelka2/strelka2_vep_tbi}
  strelka2_prepass_vcf: {type: File, outputSource: run_strelka2/strelka2_prepass_vcf}
  strelka2_vep_maf: {type: File, outputSource: run_strelka2/strelka2_vep_maf}
  manta_pass_vcf: {type: File, outputSource: run_manta/manta_pass_vcf}
  manta_prepass_vcf: {type: File, outputSource: run_manta/manta_prepass_vcf}
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
  python_vardict_interval_split:
    run: ../tools/python_vardict_interval_split.cwl
    doc: "Custom interval list generation for vardict input. Briefly, ~60M bp per interval list, 20K bp intervals, lists break on chr and N reginos only"
    in:
      wgs_bed_file: wgs_calling_interval_list
    out: [split_intervals_bed]
  
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

  gatk_filter_germline:
    run: ../tools/gatk_filter_germline_variant.cwl
    in:
      input_vcf: b_allele
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
      bed_invtl_split: python_vardict_interval_split/split_intervals_bed
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

  bedops_gen_lancet_intervals:
    run: ../tools/preprocess_lancet_intervals.cwl
    in:
      strelka2_vcf: run_strelka2/strelka2_vep_vcf
      mutect2_vcf: run_mutect2/mutect2_vep_vcf
      ref_bed: lancet_calling_interval_bed
      output_basename: output_basename
    out: [run_bed]

  gatk_intervallisttools_exome_plus:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: bedops_gen_lancet_intervals/run_bed
      reference_dict: reference_dict
      exome_flag:
        valueFrom: ${return "Y";}
      scatter_ct:
        valueFrom: ${return 50}
      bands:
        valueFrom: ${return 80000000}
    out: [output]

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
      bed_invtl_split: gatk_intervallisttools_exome_plus/output
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
      wgs_mode: cnvkit_wgs_mode
      b_allele_vcf: gatk_filter_germline/filtered_pass_vcf
      annotation_file: cnvkit_annotation_file
      output_basename: output_basename
      sex: cnvkit_sex
    out:
      [cnvkit_cnr, cnvkit_cnn_output, cnvkit_calls, cnvkit_metrics, cnvkit_gainloss, cnvkit_seg]

  run_theta2_purity:
    run: ../sub_workflows/kfdrc_run_theta2_sub_wf.cwl
    in:
      tumor_cns: run_cnvkit/cnvkit_calls
      reference_cnn: run_cnvkit/cnvkit_cnn_output
      tumor_sample_name: input_tumor_name
      normal_sample_name: input_normal_name
      paired_vcf: run_vardict/vardict_prepass_vcf
      combined_include_expression: combined_include_expression
      combined_exclude_expression: combined_exclude_expression
      min_theta2_frac: min_theta2_frac
      output_basename: output_basename
    out:
      [theta2_adjusted_cns, theta2_adjusted_seg, theta2_subclonal_results, theta2_subclonal_cns, theta2_subclone_seg]

  run_manta:
    run: ../sub_workflows/kfdrc_manta_sub_wf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      reference_dict: reference_dict
      hg38_strelka_bed: hg38_strelka_bed
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      vep_cache: vep_cache
      output_basename: output_basename
      select_vars_mode: select_vars_mode
    out:
      [manta_prepass_vcf, manta_pass_vcf]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 6
  - class: 'sbg:AWSInstanceType'
    value: c5.9xlarge
