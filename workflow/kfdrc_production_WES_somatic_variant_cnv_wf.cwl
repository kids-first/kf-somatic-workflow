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
  cfree_threads: {type: ['null', int], doc: "For ControlFreeC.  Recommend 16 max, as I/O gets saturated after that losing any advantage.", default: 16}
  vardict_min_vaf: {type: ['null', float], doc: "Min variant allele frequency for vardict to consider.  Recommend 0.05", default: 0.05}
  select_vars_mode: {type: string, doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression"}
  vardict_cpus: {type: ['null', int], default: 9}
  vardict_ram: {type: ['null', int], default: 18, doc: "In GB"}
  lancet_ram: {type: ['null', int], default: 12, doc: "Adjust in rare circumstances in which 12 GB is not enough"}
  vep_cache: {type: File, label: tar gzipped cache from ensembl/local converted cache}
  padded_capture_regions: {type: ['null', File], doc: "Recommend 100bp pad, for somatic variant"}
  unpadded_capture_regions: {type: ['null', File], doc: "Capture regions with NO padding for cnv calling"}
  cfree_chr_len: {type: File, doc: "file with chromosome lengths"}
  cfree_ploidy: {type: 'int[]', doc: "Array of ploidy possibilities for ControlFreeC to try"}
  cfree_mate_orientation_sample: {type: ['null', {type: enum, name: mate_orientation_sample, symbols: ["0", "FR", "RF", "FF"]}], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  cfree_mate_orientation_control: {type: ['null', {type: enum, name: mate_orientation_control, symbols: ["0", "FR", "RF", "FF"]}], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  b_allele: {type: ['null', File], secondaryFiles: ['.tbi'],  doc: "germline calls, needed for BAF.  GATK HC VQSR input recommended.  Tool will prefilter for germline and pass if expression given"}
  cnvkit_cnn_input: {type: ['null', File], doc: "If running using an existing .cnn, supply here"}
  mutect2_wgs_calling_interval_list: File
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
  ctrlfreec_bam_seg: {type: File, outputSource: run_controlfreec/ctrlfreec_ratio2seg}
  ctrlfreec_baf: {type: File, outputSource: run_controlfreec/ctrlfreec_baf}
  ctrlfreec_info: {type: File, outputSource: run_controlfreec/ctrlfreec_info}
  cnvkit_cnr: {type: File, outputSource: cnvkit/output_cnr}
  cnvkit_cnn_output: {type: ['null', File], outputSource: cnvkit/output_cnn}
  cnvkit_calls: {type: File, outputSource: cnvkit/output_calls}
  cnvkit_metrics: {type: File, outputSource: cnvkit/output_metrics}
  cnvkit_gainloss: {type: File, outputSource: cnvkit/output_gainloss}
  cnvkit_seg: {type: File, outputSource: cnvkit/output_seg}
  theta2_calls: {type: File, outputSource: cnvkit_import_theta2/theta2_adjusted_cns}
  theta2_seg: {type: File, outputSource: cnvkit_import_theta2/theta2_adjusted_seg}
  theta2_subclonal_results: {type: 'File[]', outputSource: [run_theta2/n3_graph, run_theta2/n2_results, run_theta2/best_results]}
  theta2_subclonal_cns: {type: 'File[]', outputSource: cnvkit_import_theta2/theta2_subclone_cns}
  theta2_subclone_seg: {type: 'File[]', outputSource: cnvkit_import_theta2/theta2_subclone_seg}
  strelka2_vep_vcf: {type: File, outputSource: vep_annot_strelka2/output_vcf}
  strelka2_vep_tbi: {type: File, outputSource: vep_annot_strelka2/output_tbi}
  strelka2_prepass_vcf: {type: File, outputSource: rename_strelka_samples/reheadered_vcf}
  strelka2_vep_maf: {type: File, outputSource: vep_annot_strelka2/output_maf}
  mutect2_vep_vcf: {type: File, outputSource: vep_annot_mutect2/output_vcf}
  mutect2_vep_tbi: {type: File, outputSource: vep_annot_mutect2/output_tbi}
  mutect2_prepass_vcf: {type: File, outputSource: filter_mutect2_vcf/filtered_vcf}
  mutect2_vep_maf: {type: File, outputSource: vep_annot_mutect2/output_maf}
  vardict_vep_somatic_only_vcf: {type: File, outputSource: run_vardict/vardict_vep_somatic_only_vcf}
  vardict_vep_somatic_only_tbi: {type: File, outputSource: run_vardict/vardict_vep_somatic_only_tbi}
  vardict_vep_somatic_only_maf: {type: File, outputSource: run_vardict/vardict_vep_somatic_only_maf}
  vardict_prepass_vcf: {type: File, outputSource: run_vardict/vardict_prepass_vcf}
  lancet_vep_vcf: {type: File, outputSource: run_lancet/lancet_vep_vcf}
  lancet_vep_tbi: {type: File, outputSource: run_lancet/lancet_vep_tbi}
  lancet_vep_maf: {type: File, outputSource: run_lancet/lancet_vep_maf}
  lancet_prepass_vcf: {type: File, outputSource: run_lancet/lancet_prepass_vcf}

steps:
  gatk_intervallisttools:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: wgs_calling_interval_list
      reference_dict: reference_dict
      exome_flag:
        valueFrom: ${return "Y";}
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
        type: ['null', int]
        default: 16
      reference: indexed_reference_fasta
    out: [bam_file]

  samtools_cram2bam_plus_calmd_normal:
    run: ../tools/samtools_cram2bam_plus_calmd.cwl
    in:
      input_reads: input_normal_aligned
      threads:
        type: ['null', int]
        default: 16
      reference: indexed_reference_fasta
    out: [bam_file]

  run_vardict:
    run: ../sub_workflows/kfdrc_vardict_WES_sub_wf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_tumor_aligned: samtools_cram2bam_plus_calmd_tumor/bam_file
      input_tumor_name: input_tumor_name
      input_normal_aligned: samtools_cram2bam_plus_calmd_normal/bam_file
      input_normal_name: input_normal_name
      output_basename: output_basename
      reference_dict: reference_dict
      bed_invtl_split: gatk_intervallisttools/output
      min_vaf: min_vaf
      select_vars_mode: select_vars_mode
      cpus: vardict_cpus
      ram: vardict_ram
      vep_cache: vep_cache
    out:
      [vardict_vep_somatic_only_vcf, vardict_vep_somatic_only_tbi, vardict_vep_somatic_only_maf, vardict_prepass_vcf]

  run_lancet:
    run: ../sub_workflows/kfdrc_lancet_WES_sub_wf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_tumor_aligned: samtools_cram2bam_plus_calmd_tumor/bam_file
      input_tumor_name: input_tumor_name
      input_normal_aligned: samtools_cram2bam_plus_calmd_normal/bam_file
      input_normal_name: input_normal_name
      output_basename: output_basename
      reference_dict: reference_dict
      bed_invtl_split: gatk_intervallisttools/output
      ram: lancet_ram
      window:
        valueFrom: ${return 600}
      padding:
        valueFrom: ${return 0}
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
      b_allele: b_allele
      chr_len: cfree_chr_len
      coeff_var: cfree_coeff_var
      contamination_adjustment: cfree_contamination_adjustment
      cfree_sex: cfree_sex
    out:
      [ctrlfreec_pval, ctrlfreec_config, ctrlfreec_pngs, ctrlfreec_bam_ratio, ctrlfreec_bam_seg, ctrlfreec_baf, ctrlfreec_info]

  run_cnvkit:
    run: ../sub_workflows/kfdrc_cnvkit_WES_sub_wf.cwl


$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4