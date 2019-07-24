cwlVersion: v1.0
class: Workflow
id: kfdrc_wes_cnv_caller_wf

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  input_tumor: { type: File, secondaryFiles: [.crai] }
  input_normal: { type: File, secondaryFiles: [.crai] }
  ref_chrs: {type: File, doc: "folder of reference chromosomes"}
  reference: {type: File, secondaryFiles: [.fai]}
  canvas_reference: {type: File, doc: "Canvas-ready kmer file"}
  chr_len: {type: File, doc: "file with chromosome lengths"}
  manifest: {type: File, doc: "Nextera manifest file"}
  b_allele_vcf: {type: File, doc: "vcf containing SNV b-alleles sites (only sites with PASS will be used)"}
  genomeSize_file: {type: File, doc: "GenomeSize.xml"}
  genome_fasta: {type: File, doc: "Genome.fa"}
  threads: int
  sample_name: string
  output_basename: string
  capture_regions: File
  exome_flag: {type: string, doc: "insert 'Y' if exome mode"}
  input_normal_name: string
  input_tumor_name: string
  vep_cache: {type: File, label: tar gzipped cache from ensembl/local converted cache}

outputs:
  ctrlfreec_cnv: {type: File, outputSource: control_free_c/output_cnv}
  ctrlfreec_cnv_bam_ratio: {type: File, outputSource: control_free_c/output_txt}
  ctrlfreec_cnv_pval: {type: File, outputSource: control_free_c_r/output_pval}
  ctrlfreec_cnv_png: {type: File, outputSource: control_free_c_viz/output_png}
  canvas_coverage_txt: {type: File, outputSource: canvas/output_txt}
  canvas_folder: {type: File, outputSource: canvas/output_folder}
  canvas_annotated_vcf: {type: File, outputSource: vep_annot_canvas/output_vcf}
  canvas_vep_tbi: {type: File, outputSource: vep_annot_canvas/output_tbi}
  canvas_maf: {type: File, outputSource: vep_annot_canvas/output_maf}
  canvas_annotated_txt: {type: File, outputSource: vep_annot_canvas/warn_txt}
  cnvkit_annotated_vcf: {type: File, outputSource: vep_annot_cnvkit/output_vcf}
  cnvkit_vep_tbi: {type: File, outputSource: vep_annot_cnvkit/output_tbi}
  cnvkit_maf: {type: File, outputSource: vep_annot_cnvkit/output_maf}
  cnvkit_calls: {type: File, outputSource: cnvkit/output_calls}
  cnvkit_scatter: {type: File, outputSource: cnvkit/output_scatter}
  cnvkit_diagram: {type: File, outputSource: cnvkit/output_diagram}
  

steps:
  samtools_tumor_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_tumor
      threads:
        valueFrom: ${return 36}
      reference: reference
    out: [bam_file]

  samtools_normal_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_normal
      threads:
        valueFrom: ${return 36}
      reference: reference
    out: [bam_file]

  gen_config:
    run: ../dev/gen_controlfreec_configfile.cwl
    in:
      tumor_bam: samtools_tumor_cram2bam/bam_file
      normal_bam: samtools_normal_cram2bam/bam_file
      capture_regions: capture_regions
      exome_flag: exome_flag
      chr_len: chr_len
      threads: threads
    out: [config_file]

  control_free_c: 
    run: ../dev/control_freec.cwl
    in: 
      tumor_bam: samtools_tumor_cram2bam/bam_file
      normal_bam: samtools_normal_cram2bam/bam_file
      capture_regions: capture_regions
      ref_chrs: ref_chrs
      chr_len: chr_len
      threads: threads
      output_basename: output_basename
      config_file: gen_config/config_file
    out: [output_txt, output_cnv]

  control_free_c_r:
    run: ../tools/control_freec_R.cwl
    in:
      cnv_bam_ratio: control_free_c/output_txt
      cnv_result: control_free_c/output_cnv
    out: [output_pval]

  control_free_c_viz:
    run: ../tools/control_freec_visualize.cwl
    in:
      output_basename: output_basename
      cnv_bam_ratio: control_free_c/output_txt
    out: [output_png]
  
  canvas:
    run: ../dev/canvas-paired-wes.cwl
    in: 
      tumor_bam: samtools_tumor_cram2bam/bam_file
      control_bam: samtools_normal_cram2bam/bam_file
      manifest: manifest
      b_allele_vcf: b_allele_vcf
      reference: canvas_reference
      genomeSize_file: genomeSize_file
      genome_fasta: genome_fasta
      filter_bed: capture_regions
      sample_name: sample_name
      output_basename: output_basename
    out: [output_vcf, output_txt, output_folder]

  cnvkit: 
    run: ../tools/cnvkit.cwl
    in:
      tumor_bam: samtools_tumor_cram2bam/bam_file
      normal_bam: samtools_normal_cram2bam/bam_file
      reference_fasta: reference
      target_regions: capture_regions
      output_reference_name: output_basename
      output_basename: output_basename
    out: [output_vcf, output_calls, output_scatter, output_diagram]

  vep_annot_canvas:
    run: ../tools/vep_vcf2maf.cwl
    in:
      input_vcf: canvas/output_vcf
      output_basename: output_basename
      tumor_id: input_tumor_name
      normal_id: input_normal_name
      tool_name:
        valueFrom: ${return "canvas_somatic"}
      reference: reference
      cache: vep_cache
    out: [output_vcf, output_tbi, output_maf, warn_txt]

  vep_annot_cnvkit:
    run: ../tools/vep_vcf2maf.cwl
    in:
      input_vcf: cnvkit/output_vcf
      output_basename: output_basename
      tumor_id: input_tumor_name
      normal_id: input_normal_name
      tool_name:
        valueFrom: ${return "cnvkit_somatic"}
      reference: reference
      cache: vep_cache
    out: [output_vcf, output_tbi, output_maf, warn_txt]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4
    