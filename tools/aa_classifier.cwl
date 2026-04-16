cwlVersion: v1.2
class: CommandLineTool
id: aa-classifier
doc: "This tools takes output from amplicon architect and classifies them"
requirements: 
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement 
  - class: DockerRequirement
    dockerPull: 'jluebeck/prepareaa:v0.1203.10'
  - class: ResourceRequirement
    ramMin: 4000
    coresMin: 8
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.data_repo)
      - entryname: setup_vars.sh
        entry: |-
          export AA_DATA_REPO=$(inputs.data_repo.path)
  - class: EnvVarRequirement
    envDef:
      AA_SRC: '/home/programs/AmpliconArchitect-master/src'
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      . ./setup_vars.sh
      && /home/programs/AmpliconClassifier-main/amplicon_classifier.py
inputs:
  data_repo: { type: Directory, doc: "Un-tarred reference obtained from https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/" }
  data_ref_version: { type: ['null', {type: enum, name: wgs_mode, symbols: ["GRCh38", "hg19", "GRCh37", "mm10", "GRCm38"]}], doc: "Genome version in data repo to use", default: "GRCh38", inputBinding: { position: 1, prefix: "--ref"} }
  cycles: { type: File, doc: "AA cycles file", inputBinding: { position: 1, prefix: "--cycles" } }
  graph: { type: File, doc: "AA graph file", inputBinding: { position: 1, prefix: "--graph"} }
outputs:
  amplicon_classification_profiles: 
    type: File
    outputBinding:
      glob: '*_amplicon_classification_profiles.tsv'
  gene_list:
    type: File
    outputBinding:
      glob: '*_gene_list.tsv'
