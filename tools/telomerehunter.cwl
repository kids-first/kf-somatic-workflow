cwlVersion: v1.2
class: CommandLineTool
id: telomerehunter
doc: |
  TelomereHunter is a tool for estimating telomere content from human
  whole-genome sequencing data. It is designed to take BAM files from a tumor and
  a matching control sample as input. However, it is also possible to run
  TelomereHunter with one input file. TelomereHunter extracts and sorts telomeric
  reads from the input sample(s). For the estimation of telomere content, GC
  biases are taken into account. Finally, the results of TelomereHunter are
  visualized in several diagrams.
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'quay.io/wtsicgp/cgp-telomerehunter:1.1.0'
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
    ramMin: $(inputs.ram * 1000)
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      1>&2 telomerehunter -o .
  - position: 10
    shellQuote: false
    prefix: "&&"
    valueFrom: >-
      find $(inputs.patient_id)/control_TelomerCnt_$(inputs.patient_id)/ -type f -maxdepth 1 | sed 'p;s/\(.*\)$(inputs.patient_id)/\1$(inputs.patient_id)_normal/' | xargs -n2 mv
  - position: 20
    shellQuote: false
    prefix: "&&"
    valueFrom: >-
      find $(inputs.patient_id)/tumor_TelomerCnt_$(inputs.patient_id)/ -type f -maxdepth 1 | sed 'p;s/\(.*\)$(inputs.patient_id)/\1$(inputs.patient_id)_tumor/' | xargs -n2 mv

inputs:
  input_bam_tumor: {type: File, inputBinding: { position: 2, prefix: "--inputBamTumor" }, secondaryFiles: [{pattern: ".bai", required: false}, {pattern: "^.bai", required: false}], doc: "Path to the indexed input BAM file of the tumor sample" }
  input_bam_control: {type: File, inputBinding: { position: 2, prefix: "--inputBamControl" }, secondaryFiles: [{pattern: ".bai", required: false}, {pattern: "^.bai", required: false}], doc: "Path to the indexed input BAM file of the control sample" } 
  patient_id: {type: 'string', inputBinding: { position: 2, prefix: "--pid" }, doc: "Sample name used in output files and diagrams" }
  banding_file: {type: 'File?', inputBinding: { position: 2, prefix: "--bandingFile" }, doc: "Path to a tab-separated file with information on chromosome banding. The first four columns of the table have to contain the chromosome name, the start and end position and the band name. The table should not have a header. If no banding file is specified, the banding information of hg19 will be used" }
  repeat_threshold: { type: 'int?', inputBinding: { position: 2, prefix: "--repeatThreshold" }, doc: "The number of repeats needed for a read to be classified as telomeric. If no repeat threshold is defined, TelomereHunter will calculate the repeat_threshold depending on the read length with the following formula: repeat_threshold = floor(read_length * 6/100)" }
  per_read_length: { type: 'boolean?', inputBinding: { position: 2, prefix: "--perReadLength" }, doc: "Repeat threshold is set per 100 bp read length. The used repeat threshold will be: floor(read_length * repeat_threshold/100) E.g. Setting -rt 8 -rl means that 8 telomere repeats are required per 100 bp read length. If the read length is 50 bp, the threshold is set to 4" }
  mapq_threshold: { type: 'int?', inputBinding: { position: 2, prefix: "--mappingQualityThreshold" }, doc: "The mapping quality needed for a read to be considered as mapped (default = 8)" }
  remove_duplicates: { type: 'boolean?', inputBinding: { position: 2, prefix: "--removeDuplicates" }, doc: "Reads marked as duplicates in the input bam file(s) are removed in the filtering step" }
  repeats: { type: 'string[]?', inputBinding: { position: 2, prefix: "--repeats" }, doc: "List of telomere repeat types to search for. Reverse complements are automatically generated and do not need to be specified! By default, TelomereHunter searches for t-, g-, c- and j-type repeats (TTAGGG TGAGGG TCAGGG TTGGGG)" }
  consecutive: { type: 'boolean?', inputBinding: { position: 2, prefix: "--consecutive" }, doc: "Search for consecutive repeats" }
  lower_gc: { type: 'int?', inputBinding: { position: 2, prefix: "--lowerGC" }, doc: "Lower limit used for GC correction of telomere content. The value must be an integer between 0 and 100 (default = 48)" }
  upper_gc: { type: 'int?', inputBinding: { position: 2, prefix: "--upperGC" }, doc: "Upper limit used for GC correction of telomere content. The value must be an integer between 0 and 100 (default = 52)" }
  no_filtering: { type: 'boolean?', inputBinding: { position: 2, prefix: "--noFiltering" }, doc: "If the filtering step of TelomereHunter has already been run previously, skip this step" }
  repeats_context: { type: 'string[]?', inputBinding: { position: 2, prefix: "--repeatsContext" }, doc: "List of telomere variant repeats for which to analyze the sequence context. Reverse complements are automatically generated and do not need to be specified! Counts for these telomere variant repeats (arbitrary and singleton context) will be added to the summary table. Default repeats: TCAGGG TGAGGG TTGGGG TTCGGG TTTGGG ATAGGG CATGGG CTAGGG GTAGGG TAAGGG)" }
  bp_context: { type: 'int?', inputBinding: { position: 2, prefix: "--bpContext" }, doc: "Number of base pairs on either side of the telomere variant repeat to investigate. Please use a number that is divisible by 6" }
  parallel: { type: 'boolean?', inputBinding: { position: 2, prefix: "--parallel" }, doc: "The filtering, sorting and estimating steps of the tumor and control sample are run in parallel. This will speed up the computation time of TelomereHunter" }
  plot_file_format:
    type:
      - 'null'
      - type: enum
        name: plot_file_format
        symbols: ["pdf", "png", "svg", "all"]
    default: "pdf"
    inputBinding:
      prefix: "--plotFileFormat"
      position: 2
    doc: |
      File format of output diagrams. Choose from pdf (default), png, svg or all (pdf, png and svg).
  plot_chr: { type: 'boolean?', inputBinding: { position: 2, prefix: "--plotChr" }, doc: "Make diagrams with telomeric reads mapping to each chromosome. If none of the options p1/p2/p3/p4/p5/p6/p7/p8 are chosen, all diagrams will be created" }
  plot_fractions: { type: 'boolean?', inputBinding: { position: 2, prefix: "--plotFractions" }, doc: "Make a diagram with telomeric reads in each fraction (intrachromosomal, subtelomeric, junction spanning, intratelomeric). If none of the options p1/p2/p3/p4/p5/p6/p7/p8 are chosen, all diagrams will be created" }
  plot_tel_content: { type: 'boolean?', inputBinding: { position: 2, prefix: "--plotTelContent" }, doc: "Make a diagram with the gc corrected telomere content in the analyzed samples. If none of the options p1/p2/p3/p4/p5/p6/p7/p8 are chosen, all diagrams will be created" }
  plot_gc: { type: 'boolean?', inputBinding: { position: 2, prefix: "--plotGC" }, doc: "Make a diagram with GC content distributions in all reads and in intratelomeric reads. If none of the options p1/p2/p3/p4/p5/p6/p7/p8 are chosen, all diagrams will be created" }
  plot_repeat_freq: { type: 'boolean?', inputBinding: { position: 2, prefix: "--plotRepeatFreq" }, doc: "Make histograms of the repeat frequencies per intratelomeric read. If none of the options p1/p2/p3/p4/p5/p6/p7/p8 are chosen, all diagrams will be created" }
  plot_tvr: { type: 'boolean?', inputBinding: { position: 2, prefix: "--plotTVR" }, doc: "Make plots for telomere variant repeats. If none of the options p1/p2/p3/p4/p5/p6/p7/p8 are chosen, all diagrams will be created" }
  plot_singleton: { type: 'boolean?', inputBinding: { position: 2, prefix: "--plotSingleton" }, doc: "Make plots for singleton telomere variant repeats. If none of the options p1/p2/p3/p4/p5/p6/p7/p8 are chosen, all diagrams will be created" }
  plot_none: { type: 'boolean?', inputBinding: { position: 2, prefix: "--plotNone" }, doc: "Do not make any diagrams. If none of the options p1/p2/p3/p4/p5/p6/p7/p8 are chosen, all diagrams will be created" }
  plot_rev_compl: { type: 'boolean?', inputBinding: { position: 2, prefix: "--plotRevCompl" }, doc: "Distinguish between forward and reverse complement telomere repeats in diagrams" }

  cpu: { type: 'int?', default: 8, doc: "Minimum reserved number of CPU cores for the task." }
  ram: { type: 'int?', default: 16, doc: "GB of RAM to allocate to the task." }

outputs:
  telomerehunter_merged_plots:
    type: File
    outputBinding:
      glob: '$(inputs.patient_id)/$(inputs.patient_id)_all_plots_merged.$(inputs.plot_file_format)'
  telomerehunter_summary:
    type: File
    outputBinding:
      glob: '$(inputs.patient_id)/$(inputs.patient_id)_summary.tsv'
  telomerehunter_control_telomere_cnt:
    type: File[]
    outputBinding:
      glob: '$(inputs.patient_id)/control_TelomerCnt_$(inputs.patient_id)/$(inputs.patient_id)*'
  telomerehunter_tumor_telomere_cnt:
    type: File[]
    outputBinding:
      glob: '$(inputs.patient_id)/tumor_TelomerCnt_$(inputs.patient_id)/$(inputs.patient_id)*'
