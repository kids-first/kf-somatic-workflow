cwlVersion: v1.2
class: CommandLineTool
id: ampliconsuite-pipeline
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'jluebeck/ampliconsuite-pipeline:v1.5.2'
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.cpu)
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.data_repo)
      - $(inputs.mosek_license_file)
  - class: EnvVarRequirement
    envDef:
      AA_SRC: '/home/programs/AmpliconArchitect-master/src'
      AC_SRC: '/home/programs/AmpliconClassifier-main'
      NCM_HOME: '/home/programs/NGSCheckMate-master'
      AA_DATA_REPO: $(runtime.outdir)/$(inputs.data_repo.basename)
      MOSEKLM_LICENSE_FILE: $(runtime.outdir)
      REF_CACHE: $(runtime.outdir)/cache/%2s/%2s/%s
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      env >&2
      $(inputs.cram_reference ? "&& /usr/bin/seq_cache_populate.pl -root cache " + inputs.cram_reference.path : "")
      && /home/programs/AmpliconSuite-pipeline-master/AmpliconSuite-pipeline.py --output_directory results
inputs:
  data_repo: { type: Directory, doc: "Un-tarred reference obtained from https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/" }
  ref_version: { type: 'string?', doc: "Reference genome version. Autodetected unless fastqs given as input." }
  mosek_license_file: { type: File, doc: "This tool uses some software that requires a license file. You can get a personal or institutional one from https://www.mosek.com/license/request/." }

  sorted_tumor_bam: { type: 'File?', inputBinding: {position: 2, prefix: "--bam"}, secondaryFiles: [{"pattern":".bai", "required": false}, {"pattern":"^.bai", "required": false}, {"pattern":".crai", "required": false}, {"pattern":"^.crai", "required": false}], doc: "Coordinate sorted tumor BAM file (aligned to an AA-supported reference.)" }
  normal_bam: { type: 'File?', inputBinding: {position: 2, prefix: "--normal_bam"}, secondaryFiles: [{"pattern":".bai", "required": false}, {"pattern":"^.bai", "required": false}, {"pattern":".crai", "required": false}, {"pattern":"^.crai", "required": false}], doc: "matched normal bam when running CNVKit" }
  cram_reference: { type: 'File?', doc: "Reference FASTA used to create CRAM inputs" }
  cnv_bed: { type: 'File?', inputBinding: {position: 2, prefix: "--cnv_bed"}, doc: "BED file (or CNVKit .cns file) of CNV changes. Fields in the bed file should be: chr start end name cngain" }
  sv_vcf: { type: 'File?', inputBinding: {position: 2, prefix: "--sv_vcf"}, doc: "Provide a VCF file of externally-called SVs to augment SVs identified by AA internally" }

  sample_name: { type: 'string?', inputBinding: {position: 2, prefix: "--sample_name"}, doc: "Sample name" }

  extra_args: { type: 'string?', inputBinding: {position: 9, shellQuote: false}, doc: "Extra args for this task" }

#  run_aa: { type: 'boolean?', doc: "Run AA after all files prepared. Default off." }
#  run_ac: { type: 'boolean?', doc: "Run AmpliconClassifier after all files prepared. Default off." }
#  no_qc: { type: 'boolean?', doc: "Skip QC on the BAM file. Do not adjust AA insert_sdevs for poor-quality insert size distribution" }
#  align_only: { type: 'boolean?', doc: "Only perform the alignment stage (do not run CNV calling and seeding)" }
#
#  download_repo: { type: 'string?', doc: "Download the selected data repo to the $AA_DATA_REPO directory and exit. '_indexed' suffix indicates BWA index is included, which is useful if performing alignment with AmpliconSuite-pipeline, but has a larger filesize." }
#
#  cngain: { type: 'float?', doc: "CN gain threshold to consider for AA seeding" }
#  cnsize_min: { type: 'int?', doc: "CN interval size (in bp) to consider for AA seeding" }
#  no_filter: { type: 'boolean?', doc: "Do not run amplified_intervals.py to remove low confidence candidate seed regions overlapping repetitive parts of the genome" }

  cpu: { type: 'int?', default: 16, inputBinding: {position: 2, prefix: "--nthreads"}, doc: 'CPUs to allocate to this task'}
  ram: { type: 'int?', default: 32, doc: 'GB of RAM to allocate to this task' }
outputs:
  log:
    type: File
    outputBinding:
      glob: 'results/*.log'
  timing_log:
    type: File
    outputBinding:
      glob: 'results/*_timing_log.txt'
  run_meta:
    type: File
    outputBinding:
      glob: 'results/*_run_metadata.json'
  sample_meta:
    type: File
    outputBinding:
      glob: 'results/*_sample_metadata.json'
  seeds:
    type: 'File?'
    outputBinding:
      glob: 'results/*_AA_CNV_SEEDS.bed'
  cnvkit_results:
    type: 'Directory?'
    outputBinding:
      glob: 'results/*_cnvkit_output'
  aa_results:
    type: 'Directory?'
    outputBinding:
      glob: 'results/*_AA_results'
  aa_summary:
    type: 'File?'
    outputBinding:
      glob: 'results/*_AA_results/*_summary.txt'
  aa_graph:
    type: 'File[]?'
    outputBinding:
      glob: 'results/*_AA_results/*_graph.txt'
  aa_cycles:
    type: 'File[]?'
    outputBinding:
      glob: 'results/*_AA_results/*_cycles.txt'
  aa_sv_png:
    type: 'File[]?'
    outputBinding:
      glob: 'results/*_AA_results/*.png'
  aa_sv_pdf:
    type: 'File[]?'
    outputBinding:
      glob: 'results/*_AA_results/*.pdf'
  ac_results:
    type: 'Directory?'
    outputBinding:
      glob: 'results/*_classification'
  ac_profiles:
    type: 'File[]?'
    outputBinding:
      glob: 'results/*_classification/*_amplicon_classification_profiles.tsv'
  ac_gene_list:
    type: 'File[]?'
    outputBinding:
      glob: 'results/*_classification*_gene_list.tsv'
