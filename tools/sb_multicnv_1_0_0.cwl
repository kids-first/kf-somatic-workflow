class: CommandLineTool
cwlVersion: v1.1
id: sb_multicnv_1_0_0
label: SB MultiCNV 1.0.0
doc: |-
  In order to generate an automated solution which could work with an unlimited number of CNV callers or tumor samples, Seven Bridges has introduced additional changes to the SB Conseca CNV workflow logic. We have created a set of tools which can analyse combined results from different callers and present the best results needed for given analysis. 

  **SB MultiCNV** consists of 3 parts: 
  - **SB MultiCNV Combine**
  - **SB MultiCNV Interpret**
  - **SB MultiCNV Segment**


  1. **SB MultiCNV Combine**

  **SB MultiCNV Combine** tool is created in order to present the results from each caller. It divides regions in which callers are in disagreement and gathers results in the form of a .CSV file. In this way, a pool of data from different callers is created from which we can derive conclusions. Two types of input files should be provided in the folder which is provided on the input:
  - Resulting CNV files obtained from different callers (in this case Control FREEC, GATK4 CNV and PURPLE were used) 
  - Information about the ploidy of the tumor sample. In the example used for this report, only PURPLE had information about ploidy as an output. If this file is not provided, the ploidy parameter is set to 2.


  2. **SB MultiCNV Interpret**

  Second tool is **SB MultiCNV Interpret **which aims to explain the results obtained from **SB MultiCNV Combine** tool. Currently, there are two ways of interpretation:
  - **Precise call** - Extracts regions where all callers provide the same status
  - **Majority call** - Extracts regions where the majority of callers agree on the status. For each call which does not fit into any category a new status is defined: *ambiguous*.

  3. **SB MultiCNV Segment**

  After determining the final call status for each region, regions which have the same status should be merged using the **SB MultiCNV Segment** tool.

  **Important**

  Input files need to be named in the following way: **caseid.tumor.caller.caller.extension** and put in the **folder**. **Do not put any '_' in the names of case_d, tumor_id and caller**. File with ploidy information and resulted file need to have the same **caseid.tumor.caller** name, but ploidy file needs to have one additional **'ploidy'** string in the filename **(case.tumor.caller.ploidy.ext)**

  Currently, MultiCNV can do the analysis of output files from the following callers:
  - **PURPLE** ('.cnv.somatic.tsv')
  - **GATK** ('.called.seg', '.modelFinal.seg')
  - **PureCN** ('.seg', '_loh.csv')
  - **CNVkit** ('.no_loh.cns', '.cns')
  - **ControlFREEC** ('.value.txt')
  - **VarSimLab** ('.CNV_results.txt')
  - **CNVnator** ('calling_results.txt')
  - **Sequenza** ('.segments.txt')
  - **Facets** ('.cncf.tsv')
  - **SCNVSim** ('.scnvsim.bed')
  - **SBG Conseca CNV** ('.conseca.tsv')
  - **Sclust** ('.allelic_states.txt')

requirements:
- class: ShellCommandRequirement
- class: LoadListingRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: images.sbgenomics.com/jelisaveta_ilic/multicnv:1.0
- class: InitialWorkDirRequirement
  listing:
  - entryname: sbg_cveto_util.py
    writable: false
    entry:
      $include: ../scripts/sbg_cveto_util.py
  - entryname: sbg_multicnv_methods.py
    writable: false
    entry:
      $include: ../scripts/sbg_multicnv_methods.py
  - entryname: sbg_multicnv.py
    writable: false
    entry:
      $include: ../scripts/sbg_multicnv.py
  - entryname: sbg_multicnv_visualisations.py
    writable: false
    entry:
      $include: ../scripts/sbg_multicnv_visualisations.py

baseCommand:
- python3
- sbg_multicnv.py

inputs:
  files:
    type: File[]
    inputBinding:
      prefix: --files
      position: 0
      shellQuote: false
    doc: |-
      Files - output files from numerous callers and files with information about the ploidy.

outputs:
  benchmark_metrics:
    type: File?
    outputBinding:
      glob: benchmark_metrics.csv
    label: Benchmark metrics
    doc: |-
      When MultiCNV is used for benchmark of callers, this outputs gives overview of the performance of callers used.
    sbg:fileTypes: CSV
  combine_heatmaps:
    type: File[]?
    outputBinding:
      glob: '*.heatmap.png'
    label: Heatmaps
    doc: |-
      Two types of heatmaps:
      1.  Comparison between each two callers - every cell shows the total length of the genome for different callers’ statuses.
      2. The comparison between copy ratio values between callers. Note: in order to make a uniform matrix, all relative copy ratios higher than 6 are set to 6. Range which represents ‘neutral’ states can be considered to be from 0.8 to 1.2.
    sbg:fileTypes: PNG
  combine_status_CN_table_overview:
    type: File?
    outputBinding:
      glob: '*.status_CN.csv'
    label: Combined regions - CN table
    doc: Resulted .CSV file which present combined regions with copy number values
    sbg:fileTypes: CSV
  combine_status_CR_table_overview:
    type: File?
    outputBinding:
      glob: '*.status_CR.csv'
    label: Combined regions - CR table
    doc: Resulted .CSV file which present combined regions with copy ratio values
    sbg:fileTypes: CSV
  combine_stats_per_file:
    type: File[]?
    outputBinding:
      glob: '*.STATS_per_file.csv'
    label: Statistics per file
    doc: |-
      Statistics about each status of genome length ('gain', 'neutral', 'loss') for each file
    sbg:fileTypes: CSV
  combine_stats_total:
    type: File[]?
    outputBinding:
      glob: '*.STATS.csv'
    label: Total statistics
    doc: |-
      Total length of genome per each combination of statuses / status ('gain', 'neutral', 'loss').
    sbg:fileTypes: CSV
  combine_status_table_overview:
    type: File?
    outputBinding:
      glob: '*.status.csv'
    label: Combine regions - status table
    doc: Resulted .CSV file which present combined regions with status values
    sbg:fileTypes: CSV
  combine_uniformed_view:
    type: File[]?
    outputBinding:
      glob: '*.uniformed.csv'
    label: Uniformed view per file
    doc: Uniformed view per file with combined regions.
    sbg:fileTypes: CSV
  html_visualisations:
    type: File[]?
    outputBinding:
      glob: '*.html'
    label: Interactive visualisations
    doc: |-
      Interactive visualisations presenting:
      1. combined calls
      2. sum of lengths
      3. total sum of lengths
    sbg:fileTypes: HTML
  interpret_majority:
    type: File?
    outputBinding:
      glob: '*.interpret_majority.csv'
    label: Result for 'majority' interpretation
    doc: |-
      Resulting .CSV file produced after 'majority' interpretation method is used to combine regions.
    sbg:fileTypes: CSV
  interpret_precise:
    type: File?
    outputBinding:
      glob: '*.interpret_precise.csv'
    label: Result for 'precise' interpretation
    doc: |-
      Resulting .CSV file produced after 'precise' interpretation method is used to combine regions.
    sbg:fileTypes: CSV
  ploidy_information:
    type: File?
    outputBinding:
      glob: '*.ploidy_overview.json'
    label: Ploidy Information
    doc: Information about the ploidy of different callers / tumors.
    sbg:fileTypes: JSON
  segment_result_final:
    type: File[]?
    outputBinding:
      glob: '*.FINAL_result.csv'
    label: Result after segmentation
    doc: |-
      Resulted file after merging the regions with the same final status excluding 'ambiguous' regions'
    sbg:fileTypes: CSV
  segment_results_with_ambiguous:
    type: File[]?
    outputBinding:
      glob: '*.seg_with_ambiguous.csv'
    label: Result after segmentation with ambiguous
    doc: Resulted file after merging the regions with the same final status
    sbg:fileTypes: CSV

$namespaces:
  sbg: https://sevenbridges.com
