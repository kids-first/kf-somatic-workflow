# main cmd to work off of:
# /home/programs/PrepareAA-master/PrepareAA.py -s BS_922YMFYK -t 8 --cnv_bed a1266520-a0a7-489e-a9c0-0596563ef178_ESTIMATED_PLOIDY_CORRECTED_CN.bed --sorted_bam a1266520-a0a7-489e-a9c0-0596563ef178.bam --downsample 10 --ref GRCh38 2> errs.log > out.log  &
cwlVersion: v1.2
class: CommandLineTool
id: prepare-aa
doc: "This tool prepares input for amplicon architect anaylsis"
requirements: 
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement 
  - class: DockerRequirement
    dockerPull: 'jluebeck/prepareaa:latest'
  - class: ResourceRequirement
    ramMin: 16000
    coresMin: $(inputs.threads)
baseCommand: [tar, -xzf]
arguments: 
  - position: 2
    shellQuote: false
    valueFrom: >-
      export AA_DATA_REPO=$PWD/$(inputs.data_untar_name)
      && touch $AA_DATA_REPO/coverage.stats
arguments: 
  - position: 3
    shellQuote: false
    valueFrom: >-
      /home/programs/PrepareAA-master/PrepareAA.py

inputs:
  data_repo: { type: File, doc: "Reference tar ball obtained from https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/", inputBinding: { position: 1} }
  data_untar_name: { type: 'string?', doc: "Name of dir created when data repo un-tarred", default: "GRCh38", inputBinding: { position: 3, prefix: "--ref"} }
  sorted_bam: { type: File, doc: "tumor bam file", secondaryFiles: [^.bai], inputBinding: { position: 3, prefix: "--sorted_bam" } }
  sample: { type: 'string', doc: "Sample name", inputBinding: { position: 3, prefix: "-s"} }
  cnv_bed: { type: File, doc: "Converted CNVkit cns-to-bed file", inputBinding: { position: 3, prefix: "--cnv_bed"} }
outputs:
  aa_cnv_seeds: 
    type: File
    outputBinding:
      glob: '*_AA_CNV_SEEDS.bed'
