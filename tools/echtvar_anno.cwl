cwlVersion: v1.2
class: CommandLineTool
id: echtvar_anno
doc: "annotate a VCF/BCF with one or more echtvar files"
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "pgc-images.sbgenomics.com/d3b-bixu/echtvar:0.1.9"
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
    ramMin: $(inputs.ram * 1000)

baseCommand: []

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      echtvar anno
  - position: 10
    prefix: "&&"
    shellQuote: false
    valueFrom: >-
      bcftools index
  - position: 19
    shellQuote: false
    valueFrom: >-
      $(inputs.output_filename)

inputs:
  input_vcf: { type: File, secondaryFiles: [{ pattern: ".tbi", required: false }, { pattern: ".csi", required: false }], inputBinding: { position: 8 }, doc: "input vcf or bcf" }
  output_filename: { type: "string?", default: "annotated.vcf.gz", inputBinding: { position: 9 }, doc: "path to bcf/vcf.gz output file (determined by extension)" }
  echtvar_zips: { type: ['null', { type: 'array', items: File, inputBinding: { prefix: "-e" }}], inputBinding: { position: 2 }, doc: "echtvar files to annotate with. can be specified many times" }
  include_expression: { type: "string?", inputBinding: { position: 2, prefix: "-i" }, doc: "expression that determines which variants to keep in output" }
  cpu: { type: 'int?', default: 8, doc: "CPUs to allocate to this task." }
  ram: { type: 'int?', default: 16, doc: "GB of RAM to allocate to this task." }

  # BCFtools Index Options
  force: { type: 'boolean?', inputBinding: { position: 12, prefix: "--force"}, doc: "overwrite index if it already exists" }
  min_shift: { type: 'int?', inputBinding: { position: 12, prefix: "--min-shift"}, doc: "set minimal interval size for CSI indices to 2^INT [14]" }
  csi: { type: 'boolean?', inputBinding: { position: 12, prefix: "--csi"}, doc: "generate CSI-format index for VCF/BCF files [default]" }
  tbi: { type: 'boolean?', inputBinding: { position: 12, prefix: "--tbi"}, doc: "generate TBI-format index for VCF files" }
  nrecords: { type: 'boolean?', inputBinding: { position: 12, prefix: "--nrecords"}, doc: "print number of records based on existing index file" }
  stats: { type: 'boolean?', inputBinding: { position: 12, prefix: "--stats"}, doc: "print per contig stats based on existing index file" }

outputs:
  annotated_vcf:
    type: File
    secondaryFiles: [{ pattern: ".tbi", required: false }, { pattern: ".csi", required: false }]
    outputBinding:
      glob: $(inputs.output_filename)
