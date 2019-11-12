cwlVersion: v1.0
class: Workflow
id: mutect2_support
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict]}
  reference_dict: File
  wgs_calling_interval_list: "File[]"
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

  exac_common_vcf: {type: File, secondaryFiles: [.tbi]}
  output_basename: string
  f1r2_counts: {type: "File[]", doc: "orientation counts from mutect2 outputs"}

outputs:
  contamination_table: {type: File, outputSource: gatk_calculate_contamination/contamination_table}
  segmentation_table: {type: File, outputSource: gatk_calculate_contamination/segmentation_table}
  f1r2_bias: {type: File, outputSource: gatk_learn_orientation_bias/f1r2_bias}
  
steps:

  gatk_learn_orientation_bias:
    run: ../tools/gatk_learnorientationbias.cwl
    label: Gatk learn bias
    in:
      input_tgz: f1r2_counts
      output_basename: output_basename
      tool_name:
        valueFrom: ${return "mutect2"}
    out: [f1r2_bias]

  gatk_get_tumor_pileup_summaries:
    label: GATK tumor pileup scatter
    run: ../tools/gatk_getpileupsummaries.cwl
    in:
      aligned_reads: input_tumor_aligned
      reference: indexed_reference_fasta
      interval_list: wgs_calling_interval_list
      exac_common_vcf: exac_common_vcf
    scatter: [interval_list]
    out: [pileup_table]

  gatk_get_normal_pileup_summaries:

    label: GATK normal pileup scatter
    run: ../tools/gatk_getpileupsummaries.cwl
    in:
      aligned_reads: input_normal_aligned
      reference: indexed_reference_fasta
      interval_list: wgs_calling_interval_list
      exac_common_vcf: exac_common_vcf
    scatter: [interval_list]
    out: [pileup_table]

  gatk_gather_tumor_pileup_summaries:
    label: GATK merge tumor pileup tables
    run: ../tools/gatk_gatherpileupsummaries.cwl
    in:
      input_tables: gatk_get_tumor_pileup_summaries/pileup_table
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name:
        valueFrom: ${return "mutect2"}
    out: [merged_table]

  gatk_gather_normal_pileup_summaries:
    label: GATK merge normal pileup tables
    run: ../tools/gatk_gatherpileupsummaries.cwl
    in:
      input_tables: gatk_get_normal_pileup_summaries/pileup_table
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name:
        valueFrom: ${return "mutect2"}
    out: [merged_table]
  
  gatk_calculate_contamination:
    run: ../tools/gatk_calculatecontamination.cwl
    in:
      tumor_pileup: gatk_gather_tumor_pileup_summaries/merged_table
      normal_pileup: gatk_gather_normal_pileup_summaries/merged_table
      output_basename: output_basename
    out: [contamination_table, segmentation_table]

$namespaces:
  sbg: https://sevenbridges.com
