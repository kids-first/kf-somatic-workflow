cwlVersion: v1.2
class: CommandLineTool
id: sentieon_tn_haplotyper2
doc: "Run Mutect2 equivalient as specified here: https://support.sentieon.com/manual/TNseq_usage/tnseq/#variant-discovery-with-matched-normal-sample"
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
    ramMin: $(inputs.ram * 1000)
  - class: DockerRequirement
    dockerPull: pgc-images.sbgenomics.com/hdchen/sentieon:202503.01
  - class: EnvVarRequirement
    envDef:
      SENTIEON_LICENSE: $(inputs.sentieon_license)
  - class: InlineJavascriptRequirement

baseCommand: [sentieon]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >
      driver
  - position: 10
    shellQuote: false
    valueFrom: >
      --algo TNhaplotyper2
  - position: 20
    shellQuote: false
    valueFrom: >
      $(inputs.output_basename).TNhaplotyper2_intermediate.vcf.gz
      --algo OrientationBias
      --tumor_sample $(inputs.tumor_sample_name)
      $(inputs.output_basename).f1r2_bias
      --algo ContaminationModel
      --tumor_sample $(inputs.tumor_sample_name)
      --normal_sample $(inputs.normal_sample_name)
      --tumor_segments $(inputs.output_basename).segments
  - position: 30
    shellQuote: false
    valueFrom: >
      $(inputs.output_basename).contamination.table 

inputs:
  # Required Arguments
  sentieon_license: { type: 'string', doc: "License server host and port" }
  indexed_reference_fasta: { type: 'File', secondaryFiles: [{ pattern: '.fai', required: true }, { pattern: '^.dict', required: true }],
    inputBinding: { position: 2, prefix: "-r"}, doc: "Reference file (FASTA)", "sbg:fileTypes": "FA, FASTA"
  }
  output_basename: { type: 'string', doc: "Output basename." }

  # Driver Arguments
  input_dedup_aligned: { type: {type: array, items: File, inputBinding: {prefix: "-i"}}, secondaryFiles: [{ pattern: '.bai', required: false}, { pattern: '^.bai', required: false }, { pattern: '.crai', required: false}, { pattern: '^.crai', required: false }], inputBinding: { position: 2, shellQuote: false }, doc: "Dedupped tumor and normal alignment files (BAM/CRAM)", "sbg:fileTypes": "BAM, CRAM" }
  cpu: { type: 'int?', default: 8, doc: "Number of CPUs to allocate to this task.",
    inputBinding: { position: 2, prefix: "-t" }
  }

  # TNhaplotyper2 Arguments
  tumor_sample_name: { type: 'string?', inputBinding: { position: 12, prefix: "--tumor_sample" }, doc: "sample name used for tumor sample in Map reads to reference stage" }
  normal_sample_name: { type: 'string?', inputBinding: { position: 12, prefix: "--normal_sample" }, doc: "sample name used for normal sample in Map reads to reference stage" }
  gnomad_vcf: { type: 'File?', inputBinding: { position: 12, prefix: "--germline_vcf" }, doc: "the location of the population germline resource, i.e. gnomAD", secondaryFiles: [{ pattern: '.tbi', required: false}, { pattern: '.idx', required: false} ] }
  # ContaminationModel
  exac_common_vcf_vcf: { type: 'File?', inputBinding: { position: 21, prefix: "--vcf" }, doc: "the location of the population germline resource, i.e. exac", secondaryFiles: [{ pattern: '.tbi', required: false}, { pattern: '.idx', required: false} ] }
  ram: { type: 'int?', default: 16, doc: "GB size of RAM to allocate to this task." }

outputs:
  tmp_vcf:
    type: File
    secondaryFiles: [{ pattern: '.tbi', required: false }]
    outputBinding:
      glob: "*.TNhaplotyper2_intermediate.vcf.gz"
  contamination_table:
    type: File
    outputBinding:
      glob: "*.contamination.table"
  contamination_segments:
    type: File
    outputBinding:
      glob: "*.segments"
  orientation_priors:
    type: File
    outputBinding:
      glob: "*.f1r2_bias"
  stats:
    type: File
    outputBinding:
      glob: "*.stats"

$namespaces:
  sbg: https://sevenbridges.com