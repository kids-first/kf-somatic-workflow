cwlVersion: v1.0
class: CommandLineTool
id: kfdrc-vep-somatic-annotate
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 24000
    coresMin: 16
  - class: DockerRequirement
    dockerPull: 'vep:r93.7'
baseCommand: [mkdir ]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.species)
      && tar -I pigz -xf $(inputs.cache.path) -C $(inputs.species)
      && perl /ensembl-vep-release-93.7/vep
      --af
      --af_1kg
      --af_esp
      --af_gnomad
      --allele_number
      --assembly $(inputs.ref_build)
      --biotype
      --buffer_size 10000
      --cache
      --cache_version $(inputs.cache_version)
      --canonical
      --ccds
      --check_existing
      --dir_cache $(inputs.species)
      --domains
      --failed 1
      --fasta $(inputs.reference.path)
      --flag_pick_allele
      --fork 16
      --format vcf
      --gene_phenotype
      --hgvs
      --input_file $(inputs.input_vcf.path)
      --no_escape
      --no_progress
      --no_stats
      --merged
      --numbers
      --offline
      --output_file $(inputs.output_basename).$(inputs.tool_name).PASS.vep.vcf
      --pick_order canonical,tsl,biotype,rank,ccds,length
      --polyphen b
      --protein
      --pubmed
      --regulatory
      --shift_hgvs 1
      --sift b
      --species $(inputs.species)
      --symbol
      --total_length
      --tsl
      --uniprot
      --variant_class
      --vcf
      --xref_refseq
      && /ensembl-vep-release-93.7/htslib/bgzip $(inputs.output_basename).$(inputs.tool_name).PASS.vep.vcf
      && /ensembl-vep-release-93.7/htslib/tabix $(inputs.output_basename).$(inputs.tool_name).PASS.vep.vcf.gz

inputs:
  reference: {type: File,  secondaryFiles: [.fai], label: Fasta genome assembly with index}
  input_vcf: {type: File, secondaryFiles: [.tbi]}
  species: {type: string?, default: "homo_sapiens"}
  output_basename: string
  tool_name: string
  cache: {type: File, label: tar gzipped cache from ensembl/local converted cache}
  cache_version: {type: int?, doc: "Version being used, should match build version", default: 93}
  ref_build: {type: string?, doc: "Genome ref build used, should line up with cache.", default: "GRCh38" }


outputs:
  output_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
  output_tbi:
    type: File  
    outputBinding:
      glob: '*.vcf.gz.tbi'
