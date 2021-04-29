cwlVersion: v1.0
class: Workflow
id: kfdrc_production_lancet
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
inputs:
  # Required
  reference_fasta: { type: File }
  reference_fai: { type: 'File?' }
  reference_dict: { type: 'File?' }
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
  vep_cache: { type: File, doc: "tar gzipped cache from ensembl/local converted cache" }
  output_basename: { type: string, doc: "String value to use as basename for outputs" }
  wgs_or_wxs: { type: { type: enum, name: sex, symbols: ["WGS", "WXS"] }, doc: "Select if this run is WGS or WXS" }

  # Optional with One Default
  lancet_ram: { type: 'int?', default: 12, doc: "Adjust in rare circumstances in which 12 GB is not enough" }
  select_vars_mode: { type: ['null', { type: enum, name: select_vars_mode, symbols: ["gatk", "grep"] }], default: "gatk", doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression" }
  vep_ref_build: { type: 'string?', default: "GRCh38", doc: "Genome ref build used, should line up with cache" }

  # Optional with Multiple Defaults (handled in choose_defaults)
  exome_flag: { type: 'string?', doc: "Whether to run in exome mode for callers. Y for WXS, N for WGS" }
  lancet_window: { type: 'int?', doc: "Window size for lancet.  Recommend 500 for WGS; 600 for exome+" }
  lancet_padding: { type: 'int?', doc: "Recommend 0 if interval file padded already, half window size if not. Recommended: 0 for WXS; 300 for WGS" }

  # WGS only Fields
  wgs_calling_interval_list: { type: 'File?', doc: "GATK intervals list-style, or bed file.  Recommend canocical chromosomes with N regions removed" }
  lancet_calling_interval_bed: { type: 'File?', doc: "For WGS, highly recommended to use CDS bed, and supplement with region calls from strelka2 & mutect2.  Can still give calling list as bed if true WGS calling desired instead of exome+" }

  # WXS only Fields
  padded_capture_regions: { type: 'File?', doc: "Recommend 100bp pad, for somatic variant" }

  # Inputs from other solo callers
  strelka2_vep_vcf: { type: 'File?', doc: "VEP VCF output from Strelka2 caller" }
  mutect2_vep_vcf: { type: 'File?', doc: "VEP VCF output from Mutect2 caller" }

outputs:
  lancet_vep_vcf: { type: File, outputSource: run_lancet/lancet_vep_vcf }
  lancet_vep_tbi: { type: File, outputSource: run_lancet/lancet_vep_tbi }
  lancet_vep_maf: { type: File, outputSource: run_lancet/lancet_vep_maf }
  lancet_prepass_vcf: { type: File, outputSource: run_lancet/lancet_prepass_vcf }

steps:
  choose_defaults:
    run: ../tools/mode_defaults.cwl
    in:
      input_mode: wgs_or_wxs
      exome_flag: exome_flag
      lancet_padding: lancet_padding
      lancet_window: lancet_window
    out: [out_exome_flag,out_cnvkit_wgs_mode,out_i_flag,out_lancet_padding,out_lancet_window,out_vardict_padding]

  prepare_reference:
    run: ../sub_workflows/prepare_reference.cwl
    in:
      input_fasta: reference_fasta
      input_fai: reference_fai
      input_dict: reference_dict
    out: [indexed_fasta,reference_dict]

  select_interval_list:
    run: ../tools/mode_selector.cwl
    in:
      input_mode: wgs_or_wxs
      wgs_input: wgs_calling_interval_list
      wxs_input: padded_capture_regions
    out: [output]

  gatk_intervallisttools:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: select_interval_list/output
      reference_dict: prepare_reference/reference_dict
      exome_flag: choose_defaults/out_exome_flag
      scatter_ct:
        valueFrom: ${return 50}
      bands:
        valueFrom: ${return 80000000}
    out: [output]

  samtools_cram2bam_plus_calmd_tumor:
    run: ../tools/samtools_cram2bam_plus_calmd.cwl
    in:
      input_reads: input_tumor_aligned
      threads:
        valueFrom: ${return 16;}
      reference: prepare_reference/indexed_fasta
    out: [bam_file]

  samtools_cram2bam_plus_calmd_normal:
    run: ../tools/samtools_cram2bam_plus_calmd.cwl
    in:
      input_reads: input_normal_aligned
      threads:
        valueFrom: ${return 16;}
      reference: prepare_reference/indexed_fasta
    out: [bam_file]

  bedops_gen_lancet_intervals:
    run: ../tools/preprocess_lancet_intervals.cwl
    in:
      strelka2_vcf: strelka2_vep_vcf
      mutect2_vcf: mutect2_vep_vcf
      ref_bed: lancet_calling_interval_bed
      output_basename: output_basename
    out: [run_bed]

  gatk_intervallisttools_exome_plus:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: bedops_gen_lancet_intervals/run_bed
      reference_dict: prepare_reference/reference_dict
      exome_flag:
        valueFrom: ${return "Y";}
      scatter_ct:
        valueFrom: ${return 50}
      bands:
        valueFrom: ${return 80000000}
    out: [output]

  select_lancet_bed_inteval:
    run: ../tools/mode_selector.cwl
    in:
      input_mode: wgs_or_wxs
      wgs_input: gatk_intervallisttools_exome_plus/output
      wxs_input: gatk_intervallisttools/output
    out: [output]

  run_lancet:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge
    run: ../sub_workflows/kfdrc_lancet_sub_wf.cwl
    in:
      indexed_reference_fasta: prepare_reference/indexed_fasta
      input_tumor_aligned: samtools_cram2bam_plus_calmd_tumor/bam_file
      input_tumor_name: input_tumor_name
      input_normal_aligned: samtools_cram2bam_plus_calmd_normal/bam_file
      input_normal_name: input_normal_name
      output_basename: output_basename
      select_vars_mode: select_vars_mode
      reference_dict: prepare_reference/reference_dict
      bed_invtl_split: select_lancet_bed_inteval/output
      ram: lancet_ram
      window: choose_defaults/out_lancet_window
      padding: choose_defaults/out_lancet_padding
      vep_cache: vep_cache
      vep_ref_build: vep_ref_build
    out:
      [lancet_vep_vcf, lancet_vep_tbi, lancet_vep_maf, lancet_prepass_vcf]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4
