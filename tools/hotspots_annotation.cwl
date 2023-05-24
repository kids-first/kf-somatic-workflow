cwlVersion: v1.2
class: CommandLineTool
id: hotspots_annotation
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.cores)
  - class: DockerRequirement
    dockerPull: 'quay.io/biocontainers/pysam:0.21.0--py310h41dec4a_1'
  - class: InitialWorkDirRequirement
    listing:
      - entryname: hotspot.py
        entry:
          $include: ../scripts/hotspot.py
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      $(inputs.disable_hotspot_annotation ? ">&2 echo 'User elected to skip hotspot annotation' && exit 0;" : ">&2 python hotspot.py")
inputs:
  # Main Arguments
  disable_hotspot_annotation: { type: 'boolean?', doc: "Disable Hotspot Annotation and skip this task." }
  input_vcf: { type: 'File', inputBinding: { position: 99, prefix: '--vcf', shellQuote: false }, secondaryFiles: [.tbi], doc: "VCF file to annotate hotspots." }
  genomic_hotspots: { type: 'File[]?', inputBinding: { position: 1, prefix: '--genomic_hotspots', shellQuote: false }, doc: "Tab-delimited BED-formatted (chrom, zero-based chromStart, non-inclusive chromEnd) file(s) containing hg38 genomic positions corresponding to hotspots" }
  protein_indels: { type: 'File[]?', inputBinding: { position: 1, prefix: '--protein_indels', shellQuote: false }, doc: "Column name-labeled, tab-delimited file(s) containing HUGO-formatted protein names and VEP-formatted positions <start_aa_pos> for SNV hotspots" }
  protein_snvs: { type: 'File[]?', inputBinding: { position: 1, prefix: '--protein_snvs', shellQuote: false }, doc: "Column name-labeled, tab-delimited file(s) containing HUGO-formatted protein names and VEP-formatted positions <start_aa_pos>-<end_aa_pos> for INDEL hotspots" }
  output_basename: { type: 'string?', inputBinding: { position: 10, prefix: '--output_basename', shellQuote: false }, doc: "String to use as basename for output file" }

  # Advanced Arguments
  protein_colname: { type: 'string?', inputBinding: { position: 2, prefix: "--protein_colname" }, doc: "Overrides the column name in the protein_hotspots file(s) that contains the protein name information" }
  position_colname: { type: 'string?', inputBinding: { position: 2, prefix: "--position_colname" }, doc: "Overrides the column name in the protein_hotspots file(s) that contains the HGVSp short information" }
  csq_field: { type: 'string?', inputBinding: { position: 2, prefix: "--csq_field" }, doc: "Overrides the name of the CSQ field that matches the values found in the protein_colname of the protein_hotspots input(s)" }
  csq_pos:
    type:
      - 'null'
      - type: enum
        name: csq_pos
        symbols: ["HGVSp","Protein_position"]
    inputBinding:
      position: 2
      prefix: "--csq_pos"
    doc: |
      Overrides the name of the CSQ field that stores protein position information. Currently capable of processing Protein_position or HGVSp.
  csq_allele: { type: 'string?', inputBinding: { position: 2, prefix: "--csq_allele" }, doc: "Overrides the name of the CSQ field that stores allele number information" }
  csq_impact: { type: 'string?', inputBinding: { position: 2, prefix: "--csq_impact" }, doc: "Overrides the name of the CSQ field that stores impact information" }
  csq_class: { type: 'string?', inputBinding: { position: 2, prefix: "--csq_class" }, doc: "Overrides the name of the CSQ field that stores variant class information" }

  # Resource Control
  ram: { type: 'int?', default: 2, doc: "GB of RAM to allocate to this task." }
  cores: { type: 'int?', default: 1, doc: "CPU cores to allocate to this task." }
outputs:
  hotspots_vcf: { type: 'File', outputBinding: { glob: '*.gz', outputEval: '$(inputs.disable_hotspot_annotation ? inputs.input_vcf : self)' }, secondaryFiles: [.tbi] }
