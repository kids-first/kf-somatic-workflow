cwlVersion: v1.0
class: Workflow
id: control-free-c-germline-wf

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  input_sample: { type: File, secondaryFiles: [.crai] }
  ref_chrs: {type: File, doc: "folder of reference chromosomes"}
  reference: {type: File, secondaryFiles: [.fai]}
  chr_len: {type: File, doc: "file with chromosome lengths"}
  threads: int
  output_basename: string
  exome_flag: {type: string, doc: "insert 'Y' if exome mode"}

outputs:
  ctrlfreec_cnv: {type: File, outputSource: control_free_c/output_cnv}
  ctrlfreec_cnv_bam_ratio: {type: File, outputSource: control_free_c/output_txt}
  ctrlfreec_cnv_pval: {type: File, outputSource: control_free_c_r/output_pval}
  ctrlfreec_cnv_png: {type: File, outputSource: control_free_c_viz/output_png}

steps:
  samtools_tumor_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_sample
      threads:
        valueFrom: ${return 36}
      reference: reference
    out: [bam_file]

  gen_config:
    run: ../dev/gen_germline_controlfreec_config.cwl
    in:
      tumor_bam: samtools_tumor_cram2bam/bam_file
      exome_flag: exome_flag
      chr_len: chr_len
      threads: threads
    out: [config_file]

  control_free_c: 
    run: ../dev/control_freec.cwl
    in: 
      tumor_bam: samtools_tumor_cram2bam/bam_file
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

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4
    