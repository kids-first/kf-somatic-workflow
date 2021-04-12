cwlVersion: v1.0
class: CommandLineTool
id: hotspots_annotation
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: ${ return inputs.ram * 1000 }
    coresMin: $(inputs.cores)
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/hotspots:0.1.0'

baseCommand: []

arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      $(inputs.disable_hotspot_annotation ? ">&2 echo 'User elected to skip hotspot annotation' && exit 0;" : ">&2 /hotspot.py")

inputs:
  disable_hotspot_annotation: { type: 'boolean?', doc: "Disable Hotspot Annotation and skip this task." }
  input_vcf: { type: 'File', inputBinding: { position: 99, prefix: '--vcf', shellQuote: false }, secondaryFiles: [.tbi], doc: "VCF file to annotate hotspots." }
  genomic_hotspots: { type: 'File[]?', inputBinding: { position: 1, prefix: '--genomic_hotspots', shellQuote: false }, doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots" }
  protein_indels: { type: 'File[]?', inputBinding: { position: 1, prefix: '--protein_indels', shellQuote: false }, doc: "Column name-labeled, tab-delimited file(s) containing HUGO-formatted protein names and VEP-formatted positions <start_aa_pos>-<end_aa_pos> for INDEL hotspots" }
  protein_snvs: { type: 'File[]?', inputBinding: { position: 1, prefix: '--protein_snvs', shellQuote: false }, doc: "Column name-labeled, tab-delimited file(s) containing HUGO-formatted protein names and VEP-formatted positions <start_aa_pos> for SNV hotspots" }
  output_basename: { type: 'string?', inputBinding: { position: 10, prefix: '--output_basename', shellQuote: false }, doc: "String to use as basename for output file" }
  ram: { type: 'int?', default: 2, doc: "GB of RAM to allocate to this task." }
  cores: { type: 'int?', default: 1, doc: "CPU cores to allocate to this task." }
outputs:
  hotspots_vcf: { type: 'File', outputBinding: { glob: '*.gz', outputEval: '$(inputs.disable_hotspot_annotation ? inputs.input_vcf : self)' }, secondaryFiles: [.tbi] }
