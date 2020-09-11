cwlVersion: v1.0
class: Workflow
id: kfdrc_combined_somatic_wgs_cnv_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict]}
  reference_fai: File
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
  threads: {type: ['null', int], doc: "For ControlFreeC.  Recommend 16 max, as I/O gets saturated after that losing any advantage.", default: 16}
  capture_regions: {type: ['null', File], doc: "If not WGS, provide this bed file"}
  chr_len: {type: File, doc: "file with chromosome lengths"}
  ploidy: {type: 'int[]', doc: "Array of ploidy possibilities for ControlFreeC to try"}
  mate_orientation_sample: {type: ['null', {type: enum, name: mate_orientation_sample, symbols: ["0", "FR", "RF", "FF"]}], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  mate_orientation_control: {type: ['null', {type: enum, name: mate_orientation_control, symbols: ["0", "FR", "RF", "FF"]}], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  b_allele: {type: ['null', File], secondaryFiles: ['.tbi'],  doc: "germline calls, needed for BAF.  GATK HC VQSR input recommended.  Tool will prefilter for germline and pass if expression given"}
  paired_vcf: {type: File, doc: "Combined somatic and germline call file. VarDict input recommended."}
  cnvkit_cnn_input: {type: ['null', File], doc: "If running using an existing .cnn, supply here"}
  wgs_mode: {type: ['null', string], doc: "for WGS mode, input Y. leave blank for hybrid mode", default: "Y"}
  annotation_file: {type: File, doc: "refFlat.txt file"}
  coeff_var: {type: ['null', float], default: 0.05, doc: "Coefficient of variation to set window size.  Default 0.05 recommended"}
  contamination_adjustment: {type: ['null', boolean], doc: "TRUE or FALSE to have ControlFreec estimate normal contam"}
  combined_include_expression: {type: ['null', string], doc: "Filter expression if vcf has non-PASS combined calls, use as-needed"}
  combined_exclude_expression: {type: ['null', string], doc: "Filter expression if vcf has non-PASS combined calls, use as-needed"}  
  min_theta2_frac: {type: ['null', float], doc: "Minimum fraction of genome with copy umber alterations.  Default is 0.05, recommend 0.01", default: 0.01}
  cfree_sex: {type: ['null', {type: enum, name: sex, symbols: ["XX", "XY"] }], doc: "If known, XX for female, XY for male"}
  cnvkit_sex: {type: ['null', string ], doc: "If known, choices are m,y,male,Male,f,x,female,Female"}
  output_basename: string

outputs:
  ctrlfreec_pval: {type: File, outputSource: rename_outputs/ctrlfreec_pval}
  ctrlfreec_config: {type: File, outputSource: rename_outputs/ctrlfreec_config}
  ctrlfreec_pngs: {type: 'File[]', outputSource: rename_outputs/ctrlfreec_pngs}
  ctrlfreec_bam_ratio: {type: File, outputSource: rename_outputs/ctrlfreec_bam_ratio}
  ctrlfreec_bam_seg: {type: File, outputSource: convert_ratio_to_seg/ctrlfreec_ratio2seg}
  ctrlfreec_baf: {type: File, outputSource: rename_outputs/ctrlfreec_baf}
  ctrlfreec_info: {type: File, outputSource: rename_outputs/ctrlfreec_info}
  cnvkit_cnr: {type: File, outputSource: cnvkit/output_cnr}
  cnvkit_cnn_output: {type: ['null', File], outputSource: cnvkit/output_cnn}
  cnvkit_calls: {type: File, outputSource: cnvkit/output_calls}
  cnvkit_metrics: {type: File, outputSource: cnvkit/output_metrics}
  cnvkit_gainloss: {type: File, outputSource: cnvkit/output_gainloss}
  cnvkit_seg: {type: File, outputSource: cnvkit/output_seg}
  theta2_calls: {type: File?, outputSource: cnvkit_import_theta2/theta2_adjusted_cns}
  theta2_seg: {type: File?, outputSource: cnvkit_import_theta2/theta2_adjusted_seg}
  theta2_subclonal_results: {type: ['null', 'File[]'], outputSource: [run_theta2/n3_graph, run_theta2/n2_results, run_theta2/best_results]}
  theta2_subclonal_cns: {type: ['null', 'File[]'], outputSource: cnvkit_import_theta2/theta2_subclone_cns}
  theta2_subclone_seg: {type: ['null', 'File[]'], outputSource: cnvkit_import_theta2/theta2_subclone_seg}

steps:
  samtools_tumor_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_tumor_aligned
      threads:
        valueFrom: ${return 16}
      reference: indexed_reference_fasta
    out: [bam_file]

  samtools_normal_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_normal_aligned
      threads:
        valueFrom: ${return 16}
      reference: indexed_reference_fasta
    out: [bam_file]

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
      snp_vcf: gatk_filter_germline/filtered_pass_vcf
    out:
      [pileup]

  controlfreec_normal_mini_pileup:
    run: ../tools/control_freec_mini_pileup.cwl
    in:
      input_reads: samtools_normal_cram2bam/bam_file
      threads:
        valueFrom: ${return 16}
      reference: indexed_reference_fasta
      snp_vcf: gatk_filter_germline/filtered_pass_vcf
    out:
      [pileup]

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

  cnvkit: 
    run: ../tools/cnvkit_batch.cwl
    in:
      input_sample: samtools_tumor_cram2bam/bam_file
      input_control: samtools_normal_cram2bam/bam_file
      reference: indexed_reference_fasta
      annotation_file: annotation_file
      output_basename: output_basename
      wgs_mode: wgs_mode
      capture_regions: capture_regions
      b_allele_vcf: gatk_filter_germline/filtered_pass_vcf
      threads: threads
      sex: cnvkit_sex
      tumor_sample_name: input_tumor_name
      cnvkit_cnn: cnvkit_cnn_input
    out: [output_cnr, output_calls, output_scatter, output_diagram, output_metrics, output_gainloss, output_seg, output_cnn]

  bcftools_filter_combined_vcf:
    run: ../tools/bcftools_filter_vcf.cwl
    in:
      input_vcf: paired_vcf
      include_expression: combined_include_expression
      exclude_expression: combined_exclude_expression
      output_basename: output_basename
    out:
      [filtered_vcf]
  
  cnvkit_export_theta2:
    run: ../tools/cnvkit_export_theta2.cwl
    in:
      tumor_cns: cnvkit/output_calls
      reference_cnn: cnvkit/output_cnn
      paired_vcf: bcftools_filter_combined_vcf/filtered_vcf
      normal_ID: input_normal_name
      tumor_ID: input_tumor_name
    out:
      [call_interval_count, call_tumor_snp, call_normal_snp]

  run_theta2:
    run: ../tools/theta2_purity.cwl
    in:
      tumor_snp: cnvkit_export_theta2/call_tumor_snp
      normal_snp: cnvkit_export_theta2/call_normal_snp
      interval_count: cnvkit_export_theta2/call_interval_count
      output_basename: output_basename
      min_frac: min_theta2_frac
    out:
      [n2_graph, n2_results, n2_withBounds, n3_graph, n3_results, n3_withBounds, best_results]
  
  cnvkit_import_theta2:
    run: ../tools/cnvkit_import_theta2.cwl
    in:
      tumor_cns: cnvkit/output_calls
      tumor_sample_name: input_tumor_name
      output_basename: output_basename
      theta2_best_results: run_theta2/best_results
      theta2_n2_results: run_theta2/n2_results
    out:
      [theta2_adjusted_cns, theta2_adjusted_seg, theta2_subclone_cns, theta2_subclone_seg]

  
  
$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4