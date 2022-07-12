cwlVersion: v1.2
class: CommandLineTool
id: amplicon-architect
doc: "This tool prepares input for amplicon architect anaylsis"
requirements: 
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement 
  - class: DockerRequirement
    dockerPull: 'jluebeck/prepareaa:latest'
  - class: ResourceRequirement
    ramMin: 16000
    coresMin: $(inputs.threads)
baseCommand: [mkdir, data_repo]
arguments: 
  - position: 1
    shellQuote: false
    valueFrom: >- 
      && tar -C data_repo -xzf
  - position: 2
    shellQuote: false
    valueFrom: >-
      && export AA_DATA_REPO=$PWD/data_repo
      && touch $AA_DATA_REPO/coverage.stats
arguments: 
  - position: 3
    shellQuote: false
    valueFrom: >-
      /home/programs/AmpliconArchitect-master/src/AmpliconArchitect.py

inputs:
  data_repo: { type: File, doc: "Reference tar ball obtained from https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/", inputBinding: { position: 1} }
  data_ref_version: { type: ['null', {type: enum, name: wgs_mode, symbols: ["GRCh38", "hg19", "GRCh37", "mm10", "GRCm38"]}], doc: "Genome version in data repo to use", default: "GRCh38", inputBinding: { position: 3, prefix: "--ref"} }
  bam: { type: File, doc: "tumor bam file", secondaryFiles: [^.bai], inputBinding: { position: 3, prefix: "--bam" } }
  downsample: { type: 'int?', doc: "Recommended for downstream anaylsis", default: 10, inputBinding: { position: 3, prefix: "--downsample" } }
  output_basename: { type: 'string', inputBinding: { position: 3, prefix: "--out"} }
  threads: { type: 'int?', doc: 'Num threads to use', default: 8 }
  aa_bed: { type: File, doc: "Bed file from prepare-aa", inputBinding: { position: 3, prefix: "--bed"} }
outputs:
  aa_cnv_seeds: 
    type: File
    outputBinding:
      glob: '*_AA_CNV_SEEDS.bed'
