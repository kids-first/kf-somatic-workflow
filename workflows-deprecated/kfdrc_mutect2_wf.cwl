cwlVersion: v1.0
class: Workflow
id: kfdrc_mutect2_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict]}
  reference_dict: File
  wgs_calling_interval_list: File
  af_only_gnomad_vcf: {type: File, secondaryFiles: ['.tbi']}
  exac_common_vcf: {type: File, secondaryFiles: ['.tbi']}
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
  exome_flag: {type: ['null', string], doc: "set to 'Y' for exome mode"}
  vep_cache: {type: File, label: tar gzipped cache from ensembl/local converted cache}
  output_basename: string
  select_vars_mode: {type: string, doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression"}

outputs:
  mutect2_filtered_stats: {type: File, outputSource: filter_mutect2_vcf/stats_table}
  mutect2_filtered_vcf: {type: File, outputSource: filter_mutect2_vcf/filtered_vcf}
  mutect2_vep_vcf: {type: File, outputSource: vep_annot_mutect2/output_vcf}
  mutect2_vep_tbi: {type: File, outputSource: vep_annot_mutect2/output_tbi}
  mutect2_vep_maf: {type: File, outputSource: vep_annot_mutect2/output_maf}
  
steps:
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
      reference: indexed_reference_fasta
      interval_list: gatk_intervallisttools/output
      af_only_gnomad_vcf: af_only_gnomad_vcf
      exome_flag: exome_flag
    scatter: [interval_list]
    out: [mutect2_vcf, f1r2_counts, mutect_stats]

  mutect2_filter_support:
    run: ../sub_workflows/kfdrc_mutect2_filter_support_subwf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      reference_dict: reference_dict
      wgs_calling_interval_list: gatk_intervallisttools/output
      input_tumor_aligned: input_tumor_aligned
      input_normal_aligned: input_normal_aligned
      exac_common_vcf: exac_common_vcf
      output_basename: output_basename
      f1r2_counts: mutect2/f1r2_counts
    out: [contamination_table, segmentation_table, f1r2_bias]

  merge_mutect2_vcf:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    run: ../tools/gatk_mergevcfs.cwl
    label: Merge mutect2 vcf
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

  gatk_selectvariants:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.2xlarge;ebs-gp2;250
    run: ../tools/gatk_selectvariants.cwl
    label: GATK Select PASS
    in:
      input_vcf: filter_mutect2_vcf/filtered_vcf
      output_basename: output_basename
      tool_name:
        valueFrom: ${return "mutect2"}
      mode: select_vars_mode
    out: [pass_vcf]

  vep_annot_mutect2:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.4xlarge;ebs-gp2;250
    run: ../tools/vep_vcf2maf.cwl
    in:
      input_vcf: gatk_selectvariants/pass_vcf
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
