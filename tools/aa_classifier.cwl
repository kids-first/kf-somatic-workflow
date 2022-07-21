cwlVersion: v1.2
class: CommandLineTool
id: aa-classifier
doc: "This tools takes output from amplicon architect and classifies them"
requirements: 
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement 
  - class: DockerRequirement
    dockerPull: 'jluebeck/prepareaa:latest'
  - class: ResourceRequirement
    ramMin: 16000
    coresMin: $(inputs.threads)
  - class: InitialWorkDirRequirement
    listing: [$(inputs.data_repo)]
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      export AA_DATA_REPO=$(inputs.data_repo.path)
      && export AA_SRC=/home/programs/AmpliconArchitect-master/src

  - position: 1
    shellQuote: false
    valueFrom: >-
      && /home/programs/AmpliconArchitect-master/src/amplicon_classifier.py
inputs:
  data_repo: { type: Directory, doc: "Un-tarred reference obtained from https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/" }
  data_ref_version: { type: ['null', {type: enum, name: wgs_mode, symbols: ["GRCh38", "hg19", "GRCh37", "mm10", "GRCm38"]}], doc: "Genome version in data repo to use", default: "GRCh38", inputBinding: { position: 1, prefix: "--ref"} }
  cycles: { type: File, doc: "AA cycles file", inputBinding: { position: 1, prefix: "--cycles" } }
  graph: { type: File, , doc: "AA graph file", inputBinding: { position: 1, prefix: "--graph"} }
  threads: { type: 'int?', doc: 'Num threads to use', default: 8 }
outputs:
  amplicon_classification_profiles: 
    type: File
    outputBinding:
      glob: '*_amplicon_classification_profiles.tsv'
  gene_list:
    type: File
    outputBinding:
      glob: '*_gene_list.tsv'


