cwlVersion: v1.2
class: CommandLineTool
id: sentieon_tn_filter
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
      --algo TNfilter
  - position: 20
    shellQuote: false
    valueFrom: >
      $(inputs.output_basename).TNhaplotyper2.filtered.vcf.gz

inputs:
  # Required Arguments
  sentieon_license: { type: 'string', doc: "License server host and port" }
  output_basename: { type: 'string', doc: "Output basename." }

  # Driver Arguments
  cpu: { type: 'int?', default: 8, doc: "Number of CPUs to allocate to this task."}
  indexed_reference_fasta: { type: 'File', secondaryFiles: [{ pattern: '.fai', required: true }, { pattern: '^.dict', required: true }],
    inputBinding: { position: 2, prefix: "-r"}, doc: "Reference file (FASTA)", "sbg:fileTypes": "FA, FASTA"
  }

  # TNfilter Arguments
  tumor_sample_name: { type: string, inputBinding: { position: 12, prefix: "--tumor_sample" }, doc: "sample name used for tumor sample in Map reads to reference stage" }
  normal_sample_name: { type: string, inputBinding: { position: 12, prefix: "--normal_sample" }, doc: "sample name used for normal sample in Map reads to reference stage" }
  tmp_vcf: { type: 'File', inputBinding: { position: 12, prefix: "-v" }, doc: "Result vcf from TN Haplotyper 2", secondaryFiles: [{ pattern: '.tbi', required: false}, { pattern: '.idx', required: false} ] }
  contamination_table: { type: File, inputBinding: { position: 12, prefix: "--contamination" } }
  tumor_segments: { type: File, inputBinding: { position: 12, prefix: "--tumor_segments" } }
  orientation_priors: { type: File, inputBinding: { position: 12, prefix: "--orientation_priors" } }
  stats: { type: File }

  ram: { type: 'int?', default: 16, doc: "GB size of RAM to allocate to this task." }

outputs:
  filtered_vcf:
    type: File
    secondaryFiles: [{ pattern: '.tbi', required: false }]
    outputBinding:
      glob: "*.TNhaplotyper2.filtered.vcf.gz"

$namespaces:
  sbg: https://sevenbridges.com