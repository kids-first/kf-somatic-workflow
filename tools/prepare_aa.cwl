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
      && export AA_SRC=/home/programs/AmpliconArchitect-master/src
      && touch $AA_DATA_REPO/coverage.stats
      ${
        if (inputs.ref_cache){
            var cmd = " && tar -xzf " + inputs.ref_cache.path
            + " && export REF_CACHE=\"$PWD/ref/cache/%2s/%2s/%s\"";
            return cmd;
        }
        else{
            return "";
        }
      }
  - position: 3
    shellQuote: false
    valueFrom: >-
      && /home/programs/PrepareAA-master/PrepareAA.py

inputs:
  data_repo: { type: File, doc: "Reference tar ball obtained from https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/", inputBinding: { position: 1} }
  data_ref_version: { type: ['null', {type: enum, name: wgs_mode, symbols: ["GRCh38", "hg19", "GRCh37", "mm10", "GRCm38"]}], doc: "Genome version in data repo to use", default: "GRCh38", inputBinding: { position: 3, prefix: "--ref"} }
  sorted_bam: { type: File, doc: "tumor bam file", secondaryFiles: [^.bai], inputBinding: { position: 3, prefix: "--sorted_bam" } }
  sample: { type: 'string', doc: "Sample name", inputBinding: { position: 3, prefix: "-s"} }
  threads: { type: 'int?', doc: 'Num threads to use', default: 8, inputBinding: { position: 3, prefix: "-t"} }
  cnv_bed: { type: File, doc: "Converted CNVkit cns-to-bed file", inputBinding: { position: 3, prefix: "--cnv_bed"} }
  ref_cache: { type: 'File?', doc: "For cram input, provide tar ball of output of running misc/seq_cache_populate.pl from samtools on the reference fasta" }
outputs:
  aa_cnv_seeds: 
    type: File
    outputBinding:
      glob: '*_AA_CNV_SEEDS.bed'
  coverage_stats:
    type: File
    outputBinding:
      glob: 'data_repo/coverage.stats'
