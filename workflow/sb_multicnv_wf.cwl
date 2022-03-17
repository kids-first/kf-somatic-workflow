class: Workflow
cwlVersion: v1.2
id: sb_multicnv_wf
label: SB MultiCNV 1.0.0 Workflow
doc: |-
  In order to generate an automated solution which could work with an unlimited number of CNV callers or tumor samples, Seven Bridges has introduced additional changes to the SB Conseca CNV workflow logic. We have created a set of tools which can analyse combined results from different callers and present the best results needed for given analysis. 

  **SB MultiCNV** consists of 3 parts: 
  - **SB MultiCNV Combine**
  - **SB MultiCNV Interpret**
  - **SB MultiCNV Segment**



  1. **SB MultiCNV Combine**

  **SB MultiCNV Combine** tool is created in order to present the results from each caller. It divides regions in which callers are in disagreement and gathers results in the form of a .CSV file. In this way, a pool of data from different callers is created from which we can derive conclusions. Two types of input files should be provided on the input:
  - Resulting CNV files obtained from different callers
  - Information about the ploidy of the tumor sample. If this file is not provided, the ploidy parameter is set to 2.


  2. **SB MultiCNV Interpret**

  Second tool is **SB MultiCNV Interpret **which aims to explain the results obtained from **SB MultiCNV Combine** tool. Currently, there are two ways of interpretation:
  - **Precise call** - Extracts regions where all callers provide the same status
  - **Majority call** - Extracts regions where the majority of callers agree on the status. For each call which does not fit into any category a new status is defined: *ambiguous*.

  3. **SB MultiCNV Segment**

  After determining the final call status for each region, regions which have the same status are merged using the **SB MultiCNV Segment** tool.

  ### **Input files**

  **SB MultiCNV** takes on input both files with information about CNVs from caller and files with the information about the ploidy value. It is important for metadata to be set up, especially field: **Case ID** and **Sample type** (which gives information which tumor sample is being processed). If there are files with the same **Case ID** and **Sample type**, on file will be overwritten. 

  Currently, SB MultiCNV can do the analysis of CNV output files from the following callers (with appropriate file extension):
  - **PURPLE** ('.cnv.somatic.tsv')
  - **GATK** ('.called.seg', '.modelFinal.seg')
  - **PureCN** ('.seg', '_loh.csv')
  - **CNVkit** ('.no_loh.cns', '.cns')
  - **ControlFREEC** ('.value.txt')
  - **VarSimLab** ('.CNV_results.txt')
  - **CNVnator** ('calling_results.txt')
  - **Sequenza** ('.segments.txt')
  - **Facets** ('.cncf.tsv')
  - **SBG Conseca CNV** ('.conseca.tsv')
  - **Sclust** ('.allelic_states.txt'). 



  Files with information about the ploidy value can be used only for specific callers:  
    
  - **PURPLE** ('.purple.purity.tsv')
  - **ControlFREEC** ('.bam_info.txt' or '.info.txt')
  - **Facets** ('.purity_ploidy.tsv')
  - **Sclust** ('.out_cn_summary.txt'). 

    

  ### **Output files**

  | Output file name | Description |
  |------------:|:------------|
  | **Benchmark metrics** | If multiCNV is used for benchmark of CNV callers, this output gives information about recall, f-score and precision of each caller |
  | **Combine regions - status table** | Resulted .CSV file which present combined regions with status values from different callers |
  | **Combine regions - CN table** | Resulted .CSV file which present combined regions with copy numbervalues from different callers |
  | **Combine regions - CR table** | Resulted .CSV file which present combined regions with copy ratio values from different callers  |
  | **Heatmaps** | Two types of heatmaps: 1.Comparison between each two callers - every cell shows the total length of the genome for different callers’ statuses. 2.The comparison between copy ratio values between callers. Note: in order to make a uniform matrix, all relative copy ratios higher than 6 are set to 6. Range which represents ‘neutral’ states can be considered to be from 0.8 to 1.2. |
  | **Interactive visualisations** | Interactive visualisations presenting: 1.combined calls  2.sum of lengths 3.total sum of lengths. Only works with **3** input data files. |
  | **Ploidy Information** | Information about the ploidy of different callers / tumors. |
  | **Result after segmentation** | **Final ** resulting file after merging the regions with the same final status excluding 'ambiguous' regions' |
  | **Result after segmentation with ambiguous** | Resulted file after merging the regions with the same final status |
  | **Result for 'majority' interpretation** | Resulting .CSV file produced after 'majority' interpretation method is used to combine regions. |
  | **Result for 'precise' interpretation** | Resulting .CSV file produced after 'precise' interpretation method is used to combine regions. |
  | **Statistics per file** | Statistics about each status of genome length ('gain', 'neutral', 'loss') for each file |
  | **Total statistics** | Total length of genome per each combination of statuses / status ('gain', 'neutral', 'loss'). |
  | **Uniformed view per file** | Uniformed view per file with combined regions. |

requirements:
- class: ScatterFeatureRequirement
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
  input_files:
    type: File[]
  output_basename:
    type: string
  override_case_id:
    type: string?
    doc: |
      Standard MultiCNV behavior will use the file metadata to obtain the
      case_id. Provide a value here to override the metadata or if metadata
      does not exist.
  override_sample_type:
    type: string?
    doc: |
      Standard MultiCNV behavior will use the file metadata to obtain the
      sample_type. Provide a value here to override the metadata or if metadata
      does not exist.

outputs:
  benchmark_metrics:
    type: File?
    outputSource: add_basename_benchmark_metrics/output
    label: Benchmark metrics
    doc: |-
      When MultiCNV is used for benchmark of callers, this outputs gives overview of the performance of callers used.
    sbg:fileTypes: CSV
  combine_heatmaps:
    type: File[]?
    outputSource: add_basename_combine_heatmaps/output
    label: Heatmaps
    doc: |-
      Two types of heatmaps:
      1.  Comparison between each two callers - every cell shows the total length of the genome for different callers’ statuses.
      2. The comparison between copy ratio values between callers. Note: in order to make a uniform matrix, all relative copy ratios higher than 6 are set to 6. Range which represents ‘neutral’ states can be considered to be from 0.8 to 1.2.
    sbg:fileTypes: PNG
  combine_stats_per_file:
    type: File[]?
    outputSource: add_basename_combine_stats_per_file/output
    label: Statistics per file
    doc: |-
      Statistics about each status of genome length ('gain', 'neutral', 'loss') for each file
    sbg:fileTypes: CSV
  combine_stats_total:
    type: File[]?
    outputSource: add_basename_combine_stats_total/output
    label: Total statistics
    doc: |-
      Total length of genome per each combination of statuses / status ('gain', 'neutral', 'loss').
    sbg:fileTypes: CSV
  combine_status_CN_table_overview:
    type: File?
    outputSource: add_basename_combine_status_CN_table_overview/output
    label: Combined regions - CN table
    doc: Resulted .CSV file which present combined regions with copy number values
    sbg:fileTypes: CSV
  combine_status_CR_table_overview:
    type: File?
    outputSource: add_basename_combine_status_CR_table_overview/output
    label: Combined regions - CR table
    doc: Resulted .CSV file which present combined regions with copy ratio values
    sbg:fileTypes: CSV
  combine_status_table_overview:
    type: File?
    outputSource: add_basename_combine_status_table_overview/output
    label: Combine regions - status table
    doc: Resulted .CSV file which present combined regions with status values
    sbg:fileTypes: CSV
  combine_uniformed_view:
    type: File[]?
    outputSource: add_basename_combine_uniformed_view/output
    label: Uniformed view per file
    doc: Uniformed view per file with combined regions.
    sbg:fileTypes: CSV
  html_visualisations:
    type: File[]?
    outputSource: add_basename_html_visualisations/output
    label: Interactive visualisations
    doc: |-
      Interactive visualisations presenting:
      1. combined calls
      2. sum of lengths
      3. total sum of lengths
    sbg:fileTypes: HTML
  interpret_majority:
    type: File?
    outputSource: add_basename_interpret_majority/output
    label: Result for 'majority' interpretation
    doc: |-
      Resulting .CSV file produced after 'majority' interpretation method is used to combine regions.
    sbg:fileTypes: CSV
  interpret_precise:
    type: File?
    outputSource: add_basename_interpret_precise/output
    label: Result for 'precise' interpretation
    doc: |-
      Resulting .CSV file produced after 'precise' interpretation method is used to combine regions.
    sbg:fileTypes: CSV
  ploidy_information:
    type: File?
    outputSource: add_basename_ploidy_information/output
    label: Ploidy Information
    doc: Information about the ploidy of different callers / tumors.
    sbg:fileTypes: JSON
  segment_result_final:
    type: File[]?
    outputSource: add_basename_segment_result_final/output
    label: Result after segmentation
    doc: |-
      Resulted file after merging the regions with the same final status excluding 'ambiguous' regions'
    sbg:fileTypes: CSV
  segment_results_with_ambiguous:
    type: File[]?
    outputSource: add_basename_segment_results_with_ambiguous/output
    label: Result after segmentation with ambiguous
    doc: Resulted file after merging the regions with the same final status
    sbg:fileTypes: CSV

steps:
  sb_multicnv_prepare_files_for_scatter:
    run: ../tools/sb_multicnv_prepare_files_for_scatter.cwl
    label: SB MultiCNV Prepare Files For Scatter
    in:
      input_files: input_files
      override_case_id: override_case_id
      override_sample_type: override_sample_type
    out: [output_files]
  sb_multicnv_rename_files:
    run: ../tools/sb_multicnv_rename_files.cwl
    label: SB MultiCNV Rename Files
    scatter: [input_file]
    scatterMethod: dotproduct
    in:
      input_file: sb_multicnv_prepare_files_for_scatter/output_files
    out: [output_file]
  sb_multicnv_1_0_0:
    run: ../tools/sb_multicnv_1_0_0.cwl
    label: SB MultiCNV 1.0.0
    in:
      files: sb_multicnv_rename_files/output_file
    out: [ploidy_information, combine_uniformed_view, combine_stats_total, combine_stats_per_file, combine_heatmaps, segment_results_with_ambiguous, combine_status_table_overview, combine_status_CR_table_overview, combine_status_CN_table_overview, interpret_precise, interpret_majority, segment_result_final, html_visualisations, benchmark_metrics]
  add_basename_benchmark_metrics:
    run: ../tools/add_basename.cwl
    when: $(inputs.input_file != null)
    in:
      input_file: sb_multicnv_1_0_0/benchmark_metrics
      output_basename: output_basename
    out: [output]
  add_basename_combine_heatmaps:
    run: ../tools/add_basename.cwl
    when: $(inputs.input_file != null)
    scatter: [input_file]
    in:
      input_file: sb_multicnv_1_0_0/combine_heatmaps
      output_basename: output_basename
    out: [output]
  add_basename_combine_status_CN_table_overview:
    run: ../tools/add_basename.cwl
    when: $(inputs.input_file != null)
    in:
      input_file: sb_multicnv_1_0_0/combine_status_CN_table_overview
      output_basename: output_basename
    out: [output]
  add_basename_combine_status_CR_table_overview:
    run: ../tools/add_basename.cwl
    when: $(inputs.input_file != null)
    in:
      input_file: sb_multicnv_1_0_0/combine_status_CR_table_overview
      output_basename: output_basename
    out: [output]
  add_basename_combine_stats_per_file:
    run: ../tools/add_basename.cwl
    when: $(inputs.input_file != null)
    scatter: [input_file]
    in:
      input_file: sb_multicnv_1_0_0/combine_stats_per_file
      output_basename: output_basename
    out: [output]
  add_basename_combine_stats_total:
    run: ../tools/add_basename.cwl
    when: $(inputs.input_file != null)
    scatter: [input_file]
    in:
      input_file: sb_multicnv_1_0_0/combine_stats_total
      output_basename: output_basename
    out: [output]
  add_basename_combine_status_table_overview:
    run: ../tools/add_basename.cwl
    when: $(inputs.input_file != null)
    in:
      input_file: sb_multicnv_1_0_0/combine_status_table_overview
      output_basename: output_basename
    out: [output]
  add_basename_combine_uniformed_view:
    run: ../tools/add_basename.cwl
    when: $(inputs.input_file != null)
    scatter: [input_file]
    in:
      input_file: sb_multicnv_1_0_0/combine_uniformed_view
      output_basename: output_basename
    out: [output]
  add_basename_html_visualisations:
    run: ../tools/add_basename.cwl
    when: $(inputs.input_file != null)
    scatter: [input_file]
    in:
      input_file: sb_multicnv_1_0_0/html_visualisations
      output_basename: output_basename
    out: [output]
  add_basename_interpret_majority:
    run: ../tools/add_basename.cwl
    when: $(inputs.input_file != null)
    in:
      input_file: sb_multicnv_1_0_0/interpret_majority
      output_basename: output_basename
    out: [output]
  add_basename_interpret_precise:
    run: ../tools/add_basename.cwl
    when: $(inputs.input_file != null)
    in:
      input_file: sb_multicnv_1_0_0/interpret_precise
      output_basename: output_basename
    out: [output]
  add_basename_ploidy_information:
    run: ../tools/add_basename.cwl
    when: $(inputs.input_file != null)
    in:
      input_file: sb_multicnv_1_0_0/ploidy_information
      output_basename: output_basename
    out: [output]
  add_basename_segment_result_final:
    run: ../tools/add_basename.cwl
    when: $(inputs.input_file != null)
    scatter: [input_file]
    in:
      input_file: sb_multicnv_1_0_0/segment_result_final
      output_basename: output_basename
    out: [output]
  add_basename_segment_results_with_ambiguous:
    run: ../tools/add_basename.cwl
    when: $(inputs.input_file != null)
    scatter: [input_file]
    in:
      input_file: sb_multicnv_1_0_0/segment_results_with_ambiguous
      output_basename: output_basename
    out: [output]

$namespaces:
  sbg: https://sevenbridges.com
