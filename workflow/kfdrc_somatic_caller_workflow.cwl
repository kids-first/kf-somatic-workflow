cwlVersion: v1.0
class: Workflow
id: kfdrc_somatic_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict]}
  reference_dict: File
  wgs_calling_interval_list: File
  input_tumor_aligned: File
  input_normal_aligned: File
  vep_cache: {type: File, label: tar gzipped cache from ensembl/local converted cache}
  output_basename: string

outputs:
  seurat_vep_vcf: {type: File, outputSource: vep_annot_seurat/output_vcf}
  seurat_vep_html: {type: File, outputSource: vep_annot_seurat/output_html}

steps:
  gatk_intervallisttools:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: wgs_calling_interval_list
    out: [output]
  
  seurat_somatic:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.4xlarge;ebs-gp2;500
    run: ../tools/seurat.cwl
    in:
      input_tumor_bam: input_tumor_aligned
      input_normal_bam: input_normal_aligned
      interval_list: gatk_intervallisttools/output
      reference: indexed_reference_fasta
    scatter: [interval_list]
    out: [seurat_vcf]
  
  merge_seurat_vcf:
    run: ../tools/gatk_mergevcfs.cwl
    in:
      input_vcfs: seurat_somatic/seurat_vcf
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name:
        valueFrom: ${ return "seurat"}
    out: [merged_vcf]
  
  vep_annot_seurat:
    run: ../tools/variant_effect_predictor.cwl
    in:
      input_vcf: merge_seurat_vcf/merged_vcf
      output_basename: output_basename
      tool_name:
        valueFrom: ${ return "seurat"}
      reference: indexed_reference_fasta
      cache: vep_cache
    out: [output_vcf, output_html, warn_txt]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 3