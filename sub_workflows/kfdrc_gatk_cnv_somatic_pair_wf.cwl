cwlVersion: v1.1
class: Workflow
id: kfdrc_gatk_cnv_somatic_pair_wf
label: KFDRC GATK CNV Somatic Pair Workflow
doc: |
  # KFDRC GATK CNV Somatic Pair Workflow

  ![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

  The Kids First Data Resource Center GATK CNV Somatic Pair Workflow is a direct
  port of the GATK best practices workflow WDL. Workflow for running the GATK CNV
  pipeline on a matched pair (normal sample optional). Supports both WGS and WES.

  ## Generalized steps:

  - prepares a genomic intervals list with PreprocessIntervals
  - collects read coverage counts across the preprocessed intervals
  - collects counts of reference versus alternate alleles with
    CollectAllelicCounts
  - annotate the preprocessed intervals
  - denoises read coverage data against the PoN with DenoiseReadCounts using
    principal component analysis
  - plots the results of standardizing and denoising copy ratios against the PoN
  - incorporates copy ratio and allelic counts data to group contiguous copy
    ratio and allelic counts segments with ModelSegments using kernel
  segmentation and Markov-chain Monte Carlo
  - calls amplification, deletion and neutral events for the segmented copy
    ratios
  - plots the results of segmentation and estimated allele-specific copy ratios

  ## Notes from GATK:

  - The intervals argument (`input_intervals` and/or `input_interval_list`) is
    required for both WGS and WES workflows and accepts formats compatible with
  the GATK -L argument. These intervals will be padded on both sides by the
  amount specified by `padding` (default 250) and split into bins of length
  specified by `bin_length` (default 1000; specify 0 to skip binning, e.g., for
  WES).  For WGS, the intervals should simply cover the autosomal chromosomes
  (sex chromosomes may be included, but care should be taken to 1) avoid creating
  panels of mixed sex, and 2) denoise case samples only with panels containing
  only individuals of the same sex as the case samples).
  - Intervals can be blacklisted from coverage collection and all downstream
    steps by using the blacklist intervals argument (`input_exclude_intervals`
  and/or `input_exclude_intervals_list`), which accepts formats compatible with
  the GATK -XL argument. This may be useful for excluding centromeric regions,
  etc. from analysis.  Alternatively, these regions may be manually filtered from
  the final callset.
  - The sites file (common_sites) should be a Picard or GATK-style interval list.
    This is a list of sites of known variation at which allelic counts will be
  collected for use in modeling minor-allele fractions.
  - To run `gatk_funcotatesegments.cwl` you must set `run_funcotatesegments` to
    `True` and provide at least `funcotator_data_sources_tgz` and
  `funcotator_ref_version`

  ### Other Resources

  dockerfiles: https://github.com/d3b-center/bixtools
requirements:
- class: InlineJavascriptRequirement
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement

inputs:
  input_aligned_reads_tumor: {type: 'File', secondaryFiles: [{pattern: ".bai", required: false},
      {pattern: "^.bai", required: false}, {pattern: ".crai", required: false}, {
        pattern: "^.crai", required: false}], doc: "File containing aligned reads\
      \ from the tumor sample. Also include associated index.", 'sbg:fileTypes': "CRAM,BAM"}
  input_aligned_reads_normal: {type: 'File?', secondaryFiles: [{pattern: ".bai", required: false},
      {pattern: "^.bai", required: false}, {pattern: ".crai", required: false}, {
        pattern: "^.crai", required: false}], doc: "File containing aligned reads\
      \ from the matched normal sample. Also include associated index.", 'sbg:fileTypes': "CRAM,BAM"}
  reference_fasta: {type: 'File', secondaryFiles: ['.fai'], doc: "Path to reference\
      \ fasta and associated fai files.", 'sbg:fileTypes': "FASTA,FA", 'sbg:suggestedValue': {
      class: File, path: 60639014357c3a53540ca7a3, name: Homo_sapiens_assembly38.fasta,
      secondaryFiles: [{class: File, path: 60639016357c3a53540ca7af, name: Homo_sapiens_assembly38.fasta.fai}]}}
  reference_dict: {type: 'File', doc: "Path to reference dict file.", 'sbg:fileTypes': "DICT",
    'sbg:suggestedValue': {class: File, path: 60639019357c3a53540ca7e7, name: Homo_sapiens_assembly38.dict}}
  input_intervals: {type: 'string[]?', doc: "Picard or GATK-style intervals in string\
      \ format that will be analyzed."}
  input_interval_list: {type: 'File?', secondaryFiles: [{pattern: ".tbi", required: false}],
    doc: "Picard or GATK-style interval list file of regions to analyze. For WGS,\
      \ this should typically only include the autosomal chromosomes.", 'sbg:fileTypes': "INTERVALS,\
      \ INTERVAL_LIST, LIST, BED, VCF, VCF.GZ"}
  input_exclude_intervals: {type: 'string[]?', doc: "Picard or GATK-style intervals\
      \ in string format that will be ignored."}
  input_exclude_interval_list: {type: 'File?', secondaryFiles: [{pattern: ".tbi",
        required: false}], doc: "Picard or GATK-style interval list file of regions\
      \ to ignore.", 'sbg:fileTypes': "INTERVALS, INTERVAL_LIST, LIST, BED, VCF, VCF.GZ"}
  bin_length: {type: 'int?', doc: "Size of bins (in bp) for coverage collection. For\
      \ WGS, use default. For WES, set this value to 0 to skip binning."}
  padding: {type: 'int?', doc: "Amount of padding (in bp) to add to both sides of\
      \ intervals."}
  collect_counts_format: {type: ['null', {type: enum, name: 'collect_counts_format',
        symbols: ["HDF5", "TSV"]}], doc: "Desired format of outputs from CollectReadCounts."}
  common_sites: {type: 'File', secondaryFiles: [{pattern: ".tbi", required: false}],
    doc: "Picard or GATK-style interval list of common sites to use for collecting\
      \ allelic counts.", 'sbg:fileTypes': "INTERVALS, INTERVAL_LIST, LIST, BED, VCF,\
      \ VCF.GZ"}
  count_panel_of_normals: {type: 'File?', doc: "Path to read-count PoN created by\
      \ the panel workflow. Significantly reduces quality of calling if not provided!",
    'sbg:fileTypes': "HDF5"}
  output_basename: {type: 'string?', doc: "String to use as the basename for the outputs\
      \ of this pipeline. If not provided, outputs will be renamed after the aligned\
      \ reads from which they were generated."}

  # CollectAllelicCounts Advanced Options
  minimum_base_quality: {type: 'int?', doc: "Minimum base quality. Base calls with\
      \ lower quality will be filtered out of pileups."}

  # DenoiseReadCounts Advanced Options
  number_of_eigensamples: {type: 'int?', doc: "Number of eigensamples to use for denoising.\
      \ If not specified or if the number of eigensamples available in the panel of\
      \ normals is smaller than this, all eigensamples will be used."}

  # ModelSegments Advanced Options
  genotyping_base_error_rate: {type: 'float?', doc: "Maximum base-error rate for genotyping\
      \ and filtering homozygous allelic counts, if available. The likelihood for\
      \ an allelic count to be generated from a homozygous site will be integrated\
      \ from zero base-error rate up to this value. Decreasing this value will increase\
      \ the number of sites assumed to be heterozygous for modeling."}
  genotyping_homozygous_log_ratio_threshold: {type: 'float?', doc: "Log-ratio threshold\
      \ for genotyping and filtering homozygous allelic counts, if available. Increasing\
      \ this value will increase the number of sites assumed to be heterozygous for\
      \ modeling."}
  kernel_approximation_dimension: {type: 'int?', doc: "Dimension of the kernel approximation.\
      \ A subsample containing this number of data points will be used to construct\
      \ the approximation for each chromosome. If the total number of data points\
      \ in a chromosome is greater than this number, then all data points in the chromosome\
      \ will be used. Time complexity scales quadratically and space complexity scales\
      \ linearly with this parameter."}
  kernel_scaling_allele_fraction: {type: 'float?', doc: "Relative scaling S of the\
      \ kernel K_AF for allele-fraction segmentation to the kernel K_CR for copy-ratio\
      \ segmentation. If multidimensional segmentation is performed, the total kernel\
      \ used will be K_CR + S * K_AF."}
  kernel_variance_allele_fraction: {type: 'float?', doc: "Variance of Gaussian kernel\
      \ for allele-fraction segmentation, if performed. If zero, a linear kernel will\
      \ be used."}
  kernel_variance_copy_ratio: {type: 'float?', doc: "Variance of Gaussian kernel for\
      \ copy-ratio segmentation, if performed. If zero, a linear kernel will be used."}
  maximum_number_of_segments_per_chromosome: {type: 'int?', doc: "Maximum number of\
      \ segments allowed per chromosome."}
  maximum_number_of_smoothing_iterations: {type: 'int?', doc: "Maximum number of iterations\
      \ allowed for segmentation smoothing."}
  minimum_total_allele_count_case: {type: 'int?', doc: "Minimum total count for filtering\
      \ allelic counts in the case sample, if available. The default value of zero\
      \ is appropriate for matched-normal mode; increase to an appropriate value for\
      \ case-only mode."}
  minimum_total_allele_count_normal: {type: 'int?', doc: "Minimum total count for\
      \ filtering allelic counts in the matched-normal sample, if available."}
  minor_allele_fraction_prior_alpha: {type: 'float?', doc: "Alpha hyperparameter for\
      \ the 4-parameter beta-distribution prior on segment minor-allele fraction.\
      \ The prior for the minor-allele fraction f in each segment is assumed to be\
      \ Beta(alpha, 1, 0, 1/2). Increasing this hyperparameter will reduce the effect\
      \ of reference bias at the expense of sensitivity."}
  number_of_burn_in_samples_allele_fraction: {type: 'int?', doc: "Number of burn-in\
      \ samples to discard for allele-fraction model."}
  number_of_burn_in_samples_copy_ratio: {type: 'int?', doc: "Number of burn-in samples\
      \ to discard for copy-ratio model."}
  number_of_changepoints_penalty_factor: {type: 'int?', doc: "Factor A for the penalty\
      \ on the number of changepoints per chromosome for segmentation. Adds a penalty\
      \ of the form A * C * [1 + log (N / C)], where C is the number of changepoints\
      \ in the chromosome, to the cost function for each chromosome. Must be non-negative."}
  number_of_samples_allele_fraction: {type: 'int?', doc: "Total number of MCMC samples\
      \ for allele-fraction model."}
  number_of_samples_copy_ratio: {type: 'int?', doc: "Total number of MCMC samples\
      \ for copy-ratio model."}
  number_of_smoothing_iterations_per_fit: {type: 'int?', doc: "Number of segmentation-smoothing\
      \ iterations per MCMC model refit. (Increasing this will decrease runtime, but\
      \ the final number of segments may be higher. Setting this to 0 will completely\
      \ disable model refitting between iterations.)"}
  smoothing_credible_interval_threshold_allele_fraction: {type: 'float?', doc: "Number\
      \ of 10% equal-tailed credible-interval widths to use for allele-fraction segmentation\
      \ smoothing."}
  smoothing_credible_interval_threshold_copy_ratio: {type: 'float?', doc: "Number\
      \ of 10% equal-tailed credible-interval widths to use for copy-ratio segmentation\
      \ smoothing."}
  window_size: {type: 'int[]?', doc: "Window sizes to use for calculating local changepoint\
      \ costs. For each window size, the cost for each data point to be a changepoint\
      \ will be calculated assuming that the point demarcates two adjacent segments\
      \ of that size. Including small (large) window sizes will increase sensitivity\
      \ to small (large) events. Duplicate values will be ignored."}

  # CallCopyRatioSegments Advanced Options
  calling_copy_ratio_z_score_threshold: {type: 'float?', doc: "Threshold on z-score\
      \ of non-log2 copy ratio used for calling segments."}
  neutral_segment_copy_ratio_lower_bound: {type: 'float?', doc: "Lower bound on non-log2\
      \ copy ratio used for determining copy-neutral segments."}
  neutral_segment_copy_ratio_upper_bound: {type: 'float?', doc: "Upper bound on non-log2\
      \ copy ratio used for determining copy-neutral segments."}
  outlier_neutral_segment_copy_ratio_z_score_threshold: {type: 'float?', doc: "Threshold\
      \ on z-score of non-log2 copy ratio used for determining outlier copy-neutral\
      \ segments. If non-log2 copy ratio z-score is above this threshold for a copy-neutral\
      \ segment, it is considered an outlier and not used in the calculation of the\
      \ length-weighted mean and standard deviation used for calling."}

  # Plotting Advanced Options
  minimum_contig_length: {type: 'int?', doc: "Threshold length (in bp) for contigs\
      \ to be plotted. Contigs with lengths less than this threshold will not be plotted.\
      \ This can be used to filter out mitochondrial contigs, unlocalized contigs,\
      \ etc."}

  # Funcotator Advanced Options
  run_funcotatesegments: {type: 'boolean?', doc: "If true, run Funcotator on the called\
      \ copy-ratio segments. This will generate both a simple TSV and a gene list."}
  funcotator_annotation_default: {type: 'string[]?', doc: "Annotations to include\
      \ in all annotated variants if the annotation is not specified in the data sources\
      \ (in the format :). This will add the specified annotation to every annotated\
      \ variant if it is not already present."}
  funcotator_annotation_override: {type: 'string[]?', doc: "Override values for annotations\
      \ (in the format :). Replaces existing annotations of the given name with given\
      \ values."}
  funcotator_data_sources_tgz: {type: 'File?', doc: "Path to tar.gz containing the\
      \ data sources for Funcotator to create annotations.", 'sbg:fileTypes': "TAR,\
      \ TAR.GZ, TGZ", 'sbg:suggestedValue': {class: File, path: 60e5f8636a504e4e0c6408d8,
      name: funcotator_dataSources.v1.6.20190124s.tar.gz}}
  funcotator_exclude_field: {type: 'string[]?', doc: "Fields that should not be rendered\
      \ in the final output. Only exact name matches will be excluded."}
  funcotator_ref_version: {type: ['null', {type: enum, name: 'funcotator_ref_version',
        symbols: ["hg19", "hg38", "b37"]}], default: "hg38", doc: "The version of\
      \ the Human Genome reference to use (e.g. hg19, hg38, etc.). This will correspond\
      \ to a sub-folder of each data source corresponding to that data source for\
      \ the given reference."}
  funcotator_transcript_list: {type: 'string[]?', doc: "A set of transcript IDs to\
      \ use for annotation to override selected transcript."}
  funcotator_transcript_list_file: {type: 'File?', doc: "File to use as a list of\
      \ transcripts (one transcript ID per line, version numbers are ignored)"}
  funcotator_transcript_selection_mode: {type: ['null', {type: enum, name: 'funcotator_transcript_selection_mode',
        symbols: ["ALL", "BEST_EFFECT", "CANONICAL"]}], doc: "Method of detailed transcript\
      \ selection. This will select the transcript for detailed annotation (CANONICAL,\
      \ ALL, or BEST_EFFECT)."}
  funcotator_minimum_segment_size: {type: 'int?', doc: "The minimum number of bases\
      \ for a variant to be annotated as a segment. Recommended to be changed only\
      \ for use with FuncotateSegments. If you encounter 'Variant context does not\
      \ represent a copy number segment' error, set this value lower than the length\
      \ of the failed segment."}

  # Resource Control
  allelic_counts_max_mem: {type: 'int?', default: 100, doc: "Max memory in GB to allocate\
      \ to CollectAllelicCounts"}
  model_segments_max_mem: {type: 'int?', default: 30, doc: "Max memory in GB to allocate\
      \ to ModelSegments"}


outputs:
#  preprocessed_intervals: {type: File, outputSource: preprocess_intervals/output,
#    doc: "Preprocessed Picard interval-list file used to collect read counts."}

  tumor_file_archive: {type: File, outputSource: tar_outputs_tumor/output, doc: "Tar\
      \ archive containing all files generated from the tumor sample aligned reads."}
  modeled_segments_tumor: {type: File, outputSource: model_segments_tumor/modeled_segments,
    doc: "modelFinal.seg files. These are tab-separated values (TSV) files with a\
      \ SAM-style header containing a read group sample name, a sequence dictionary,\
      \ a row specifying the column headers contained in ModeledSegmentCollection.ModeledSegmentTableColumn,\
      \ and the corresponding entry rows."}
  modeled_segments_tumor_plot: {type: File, outputSource: plot_modeled_segments_tumor/output,
    doc: "Modeled-segments-plot file. This shows the input denoised copy ratios and/or\
      \ alternate-allele fractions as points, as well as box plots for the available\
      \ posteriors in each segment. The colors of the points alternate with the segmentation.\
      \ Copy ratios are only plotted up to the maximum value specified by the argument\
      \ maximum-copy-ratio. Point sizes can be specified by the arguments point-size-copy-ratio\
      \ and point-size-allele-fraction."}
  called_copy_ratio_segments_tumor: {type: File, outputSource: call_copy_ratio_segments_tumor/called_segments,
    doc: "Called copy-ratio-segments file. This is a tab-separated values (TSV) file\
      \ with a SAM-style header containing a read group sample name, a sequence dictionary,\
      \ a row specifying the column headers contained in CalledCopyRatioSegmentCollection.CalledCopyRatioSegmentTableColumn,\
      \ and the corresponding entry rows."}
  denoised_tumor_plot: {type: File, outputSource: plot_denoised_copy_ratios_tumor/denoised_plot,
    doc: "Denoised-plot file that covers the entire range of the copy ratios"}

  normal_file_archive: {type: 'File?', outputSource: tar_outputs_normal/output, doc: "Tar\
      \ archive containing all files generated from the normal sample aligned reads."}
  modeled_segments_normal: {type: 'File?', outputSource: model_segments_normal/modeled_segments,
    doc: "modelFinal.seg files. These are tab-separated values (TSV) files with a\
      \ SAM-style header containing a read group sample name, a sequence dictionary,\
      \ a row specifying the column headers contained in ModeledSegmentCollection.ModeledSegmentTableColumn,\
      \ and the corresponding entry rows."}
  modeled_segments_normal_plot: {type: 'File?', outputSource: plot_modeled_segments_normal/output,
    doc: "Modeled-segments-plot file. This shows the input denoised copy ratios and/or\
      \ alternate-allele fractions as points, as well as box plots for the available\
      \ posteriors in each segment. The colors of the points alternate with the segmentation.\
      \ Copy ratios are only plotted up to the maximum value specified by the argument\
      \ maximum-copy-ratio. Point sizes can be specified by the arguments point-size-copy-ratio\
      \ and point-size-allele-fraction."}
  called_copy_ratio_segments_normal: {type: 'File?', outputSource: call_copy_ratio_segments_normal/called_segments,
    doc: "Called copy-ratio-segments file. This is a tab-separated values (TSV) file\
      \ with a SAM-style header containing a read group sample name, a sequence dictionary,\
      \ a row specifying the column headers contained in CalledCopyRatioSegmentCollection.CalledCopyRatioSegmentTableColumn,\
      \ and the corresponding entry rows."}
  denoised_normal_plot: {type: 'File?', outputSource: plot_denoised_copy_ratios_normal/denoised_plot,
    doc: "Denoised-plot file that covers the entire range of the copy ratios"}

  funcotated_called_file_tumor: {type: 'File?', outputSource: funcotate_segments/funcotated_seg_simple_tsv,
    doc: "TSV where each row is a segment and the annotations are the covered genes\
      \ and which genes+exon is overlapped by the segment breakpoints."}
  funcotated_called_gene_list_file_tumor: {type: 'File?', outputSource: funcotate_segments/funcotated_gene_list_tsv,
    doc: "TSV where each row is a gene and the annotations are the covered genes and\
      \ which genes+exon is overlapped by the segment breakpoints."}

steps:
  preprocess_intervals:
    run: ../tools/gatk_preprocessintervals.cwl
    in:
      bin_length: bin_length
      input_intervals: input_intervals
      input_interval_list: input_interval_list
      input_exclude_intervals: input_exclude_intervals
      input_exclude_interval_list: input_exclude_interval_list
      interval_merging_rule: {default: "OVERLAPPING_ONLY"}
      padding: padding
      reference: reference_fasta
      sequence_dictionary: reference_dict
    out: [output]

  collect_read_counts_tumor:
    run: ../tools/gatk_collectreadcounts.cwl
    in:
      input_interval_list: preprocess_intervals/output
      input_aligned_reads: input_aligned_reads_tumor
      interval_merging_rule: {default: "OVERLAPPING_ONLY"}
      output_format: collect_counts_format
      reference: reference_fasta
      sequence_dictionary: reference_dict
      output_prefix:
        source: output_basename
        valueFrom: |
          $(self ? self+".tumor.gatk_cnv" : null)
    out: [output]

  collect_read_counts_normal:
    run: ../tools/gatk_collectreadcounts.cwl
    in:
      input_interval_list: preprocess_intervals/output
      input_aligned_reads: input_aligned_reads_normal
      interval_merging_rule: {default: "OVERLAPPING_ONLY"}
      output_format: collect_counts_format
      reference: reference_fasta
      sequence_dictionary: reference_dict
      output_prefix:
        source: output_basename
        valueFrom: |
          $(self ? self+".normal.gatk_cnv" : null)
    out: [output]

  collect_allelic_counts_tumor:
    run: ../tools/gatk_collectalleliccounts.cwl
    in:
      input_aligned_reads: input_aligned_reads_tumor
      input_interval_list: common_sites
      minimum_base_quality: minimum_base_quality
      reference: reference_fasta
      sequence_dictionary: reference_dict
      max_memory: allelic_counts_max_mem
      output_prefix:
        source: output_basename
        valueFrom: |
          $(self ? self+".tumor.gatk_cnv" : null)
    out: [output]

  collect_allelic_counts_normal:
    run: ../tools/gatk_collectalleliccounts.cwl
    in:
      input_aligned_reads: input_aligned_reads_normal
      input_interval_list: common_sites
      minimum_base_quality: minimum_base_quality
      reference: reference_fasta
      sequence_dictionary: reference_dict
      max_memory: allelic_counts_max_mem
      output_prefix:
        source: output_basename
        valueFrom: |
          $(self ? self+".normal.gatk_cnv" : null)
    out: [output]

  denoise_read_counts_tumor:
    run: ../tools/gatk_denoisereadcounts.cwl
    in:
      count_panel_of_normals: count_panel_of_normals
      read_counts: collect_read_counts_tumor/output
      number_of_eigensamples: number_of_eigensamples
      output_prefix:
        source: output_basename
        valueFrom: |
          $(self ? self+".tumor.gatk_cnv" : null)
    out: [denoised_copy_ratios, standardized_copy_ratios]

  denoise_read_counts_normal:
    run: ../tools/gatk_denoisereadcounts.cwl
    in:
      count_panel_of_normals: count_panel_of_normals
      read_counts: collect_read_counts_normal/output
      number_of_eigensamples: number_of_eigensamples
      output_prefix:
        source: output_basename
        valueFrom: |
          $(self ? self+".normal.gatk_cnv" : null)
    out: [denoised_copy_ratios, standardized_copy_ratios]

  model_segments_tumor:
    run: ../tools/gatk_modelsegments.cwl
    in:
      allelic_counts: collect_allelic_counts_tumor/output
      denoised_copy_ratios: denoise_read_counts_tumor/denoised_copy_ratios
      genotyping_base_error_rate: genotyping_base_error_rate
      genotyping_homozygous_log_ratio_threshold: genotyping_homozygous_log_ratio_threshold
      kernel_approximation_dimension: kernel_approximation_dimension
      kernel_scaling_allele_fraction: kernel_scaling_allele_fraction
      kernel_variance_allele_fraction: kernel_variance_allele_fraction
      kernel_variance_copy_ratio: kernel_variance_copy_ratio
      maximum_number_of_segments_per_chromosome: maximum_number_of_segments_per_chromosome
      maximum_number_of_smoothing_iterations: maximum_number_of_smoothing_iterations
      minimum_total_allele_count_case: minimum_total_allele_count_case
      minimum_total_allele_count_normal: minimum_total_allele_count_normal
      minor_allele_fraction_prior_alpha: minor_allele_fraction_prior_alpha
      normal_allelic_counts: collect_allelic_counts_normal/output
      number_of_burn_in_samples_allele_fraction: number_of_burn_in_samples_allele_fraction
      number_of_burn_in_samples_copy_ratio: number_of_burn_in_samples_copy_ratio
      number_of_changepoints_penalty_factor: number_of_changepoints_penalty_factor
      number_of_samples_allele_fraction: number_of_samples_allele_fraction
      number_of_samples_copy_ratio: number_of_samples_copy_ratio
      number_of_smoothing_iterations_per_fit: number_of_smoothing_iterations_per_fit
      smoothing_credible_interval_threshold_allele_fraction: smoothing_credible_interval_threshold_allele_fraction
      smoothing_credible_interval_threshold_copy_ratio: smoothing_credible_interval_threshold_copy_ratio
      window_size: window_size
      max_memory: model_segments_max_mem
      output_prefix:
        source: output_basename
        valueFrom: |
          $(self ? self+".tumor.gatk_cnv" : null)
    out: [allele_fraction_legacy_segments, allele_fraction_parameters, allele_fraction_parameters_begin,
      copy_ratio_legacy_segments, copy_ratio_only_segments, copy_ratio_parameters,
      copy_ratio_parameters_begin, het_allelic_counts, modeled_segments, modeled_segments_begin,
      normal_het_allelic_counts]

  model_segments_normal:
    run: ../tools/gatk_modelsegments.cwl
    in:
      allelic_counts: collect_allelic_counts_normal/output
      denoised_copy_ratios: denoise_read_counts_normal/denoised_copy_ratios
      genotyping_base_error_rate: genotyping_base_error_rate
      genotyping_homozygous_log_ratio_threshold: genotyping_homozygous_log_ratio_threshold
      kernel_approximation_dimension: kernel_approximation_dimension
      kernel_scaling_allele_fraction: kernel_scaling_allele_fraction
      kernel_variance_allele_fraction: kernel_variance_allele_fraction
      kernel_variance_copy_ratio: kernel_variance_copy_ratio
      maximum_number_of_segments_per_chromosome: maximum_number_of_segments_per_chromosome
      maximum_number_of_smoothing_iterations: maximum_number_of_smoothing_iterations
      minimum_total_allele_count_case: minimum_total_allele_count_normal
      minor_allele_fraction_prior_alpha: minor_allele_fraction_prior_alpha
      normal_allelic_counts: collect_allelic_counts_normal/output
      number_of_burn_in_samples_allele_fraction: number_of_burn_in_samples_allele_fraction
      number_of_burn_in_samples_copy_ratio: number_of_burn_in_samples_copy_ratio
      number_of_changepoints_penalty_factor: number_of_changepoints_penalty_factor
      number_of_samples_allele_fraction: number_of_samples_allele_fraction
      number_of_samples_copy_ratio: number_of_samples_copy_ratio
      number_of_smoothing_iterations_per_fit: number_of_smoothing_iterations_per_fit
      smoothing_credible_interval_threshold_allele_fraction: smoothing_credible_interval_threshold_allele_fraction
      smoothing_credible_interval_threshold_copy_ratio: smoothing_credible_interval_threshold_copy_ratio
      window_size: window_size
      max_memory: model_segments_max_mem
      output_prefix:
        source: output_basename
        valueFrom: |
          $(self ? self+".normal.gatk_cnv" : null)
    out: [allele_fraction_legacy_segments, allele_fraction_parameters, allele_fraction_parameters_begin,
      copy_ratio_legacy_segments, copy_ratio_only_segments, copy_ratio_parameters,
      copy_ratio_parameters_begin, het_allelic_counts, modeled_segments, modeled_segments_begin,
      normal_het_allelic_counts]

  call_copy_ratio_segments_tumor:
    run: ../tools/gatk_callcopyratiosegments.cwl
    in:
      calling_copy_ratio_z_score_threshold: calling_copy_ratio_z_score_threshold
      copy_ratio_segments: model_segments_tumor/copy_ratio_only_segments
      neutral_segment_copy_ratio_lower_bound: neutral_segment_copy_ratio_lower_bound
      neutral_segment_copy_ratio_upper_bound: neutral_segment_copy_ratio_upper_bound
      outlier_neutral_segment_copy_ratio_z_score_threshold: outlier_neutral_segment_copy_ratio_z_score_threshold
      output_prefix:
        source: output_basename
        valueFrom: |
          $(self ? self+".tumor.gatk_cnv" : null)
    out: [called_legacy_segments, called_segments]

  call_copy_ratio_segments_normal:
    run: ../tools/gatk_callcopyratiosegments.cwl
    in:
      calling_copy_ratio_z_score_threshold: calling_copy_ratio_z_score_threshold
      copy_ratio_segments: model_segments_normal/copy_ratio_only_segments
      neutral_segment_copy_ratio_lower_bound: neutral_segment_copy_ratio_lower_bound
      neutral_segment_copy_ratio_upper_bound: neutral_segment_copy_ratio_upper_bound
      outlier_neutral_segment_copy_ratio_z_score_threshold: outlier_neutral_segment_copy_ratio_z_score_threshold
      output_prefix:
        source: output_basename
        valueFrom: |
          $(self ? self+".normal.gatk_cnv" : null)
    out: [called_legacy_segments, called_segments]

  plot_denoised_copy_ratios_tumor:
    run: ../tools/gatk_plotdenoisedcopyratios.cwl
    in:
      denoised_copy_ratios: denoise_read_counts_tumor/denoised_copy_ratios
      minimum_contig_length: minimum_contig_length
      sequence_dictionary: reference_dict
      standardized_copy_ratios: denoise_read_counts_tumor/standardized_copy_ratios
      output_prefix:
        source: output_basename
        valueFrom: |
          $(self ? self+".tumor.gatk_cnv" : null)
    out: [denoised_plot, delta_mad, denoised_mad, scaled_delta_mad, standardized_mad]

  plot_denoised_copy_ratios_normal:
    run: ../tools/gatk_plotdenoisedcopyratios.cwl
    in:
      denoised_copy_ratios: denoise_read_counts_normal/denoised_copy_ratios
      minimum_contig_length: minimum_contig_length
      sequence_dictionary: reference_dict
      standardized_copy_ratios: denoise_read_counts_normal/standardized_copy_ratios
      output_prefix:
        source: output_basename
        valueFrom: |
          $(self ? self+".normal.gatk_cnv" : null)
    out: [denoised_plot, delta_mad, denoised_mad, scaled_delta_mad, standardized_mad]

  plot_modeled_segments_tumor:
    run: ../tools/gatk_plotmodeledsegments.cwl
    in:
      allelic_counts: model_segments_tumor/het_allelic_counts
      denoised_copy_ratios: denoise_read_counts_tumor/denoised_copy_ratios
      minimum_contig_length: minimum_contig_length
      segments: model_segments_tumor/modeled_segments
      sequence_dictionary: reference_dict
      output_prefix:
        source: output_basename
        valueFrom: |
          $(self ? self+".tumor.gatk_cnv" : null)
    out: [output]

  plot_modeled_segments_normal:
    run: ../tools/gatk_plotmodeledsegments.cwl
    in:
      allelic_counts: model_segments_normal/het_allelic_counts
      denoised_copy_ratios: denoise_read_counts_normal/denoised_copy_ratios
      minimum_contig_length: minimum_contig_length
      segments: model_segments_normal/modeled_segments
      sequence_dictionary: reference_dict
      output_prefix:
        source: output_basename
        valueFrom: |
          $(self ? self+".normal.gatk_cnv" : null)
    out: [output]

  awk_min_seg_length_tumor:
    run: ../tools/awk_min_seg_length.cwl
    in:
      input_file: call_copy_ratio_segments_tumor/called_segments
      default_min_len: funcotator_minimum_segment_size
    out: [output]

  funcotate_segments:
    run: ../tools/gatk_funcotatesegments.cwl
    in:
      run_funcotatesegments: run_funcotatesegments
      annotation_default: funcotator_annotation_default
      annotation_override: funcotator_annotation_override
      data_sources_tgz: funcotator_data_sources_tgz
      exclude_field: funcotator_exclude_field
      ref_version: {source: funcotator_ref_version, default: "hg38"}
      reference: reference_fasta
      segments: call_copy_ratio_segments_tumor/called_segments
      sequence_dictionary: reference_dict
      transcript_list: funcotator_transcript_list
      transcript_list_file: funcotator_transcript_list_file
      transcript_selection_mode: funcotator_transcript_selection_mode
      minimum_segment_size:
        source: awk_min_seg_length_tumor/output
        valueFrom: |
          $(self - 1)
      output_prefix:
        source: output_basename
        valueFrom: |
          $(self ? self+".gatk_cnv" : null)
    out: [funcotated_seg_simple_tsv, funcotated_gene_list_tsv]

  tar_outputs_tumor:
    run: ../tools/tar.cwl
    in:
      output_file:
        source: model_segments_tumor/modeled_segments
        valueFrom: |
          $(self ? self.nameroot.split('.')[0]+".tumor.gatk_cnv.tar.gz" : "tumor_outs.gatk_cnv.tar.gz")
      input_files:
        source:
        - collect_read_counts_tumor/output
        - collect_allelic_counts_tumor/output
        - denoise_read_counts_tumor/denoised_copy_ratios
        - denoise_read_counts_tumor/standardized_copy_ratios
        - model_segments_tumor/allele_fraction_legacy_segments
        - model_segments_tumor/allele_fraction_parameters
        - model_segments_tumor/allele_fraction_parameters_begin
        - model_segments_tumor/copy_ratio_legacy_segments
        - model_segments_tumor/copy_ratio_only_segments
        - model_segments_tumor/copy_ratio_parameters
        - model_segments_tumor/copy_ratio_parameters_begin
        - model_segments_tumor/het_allelic_counts
        - model_segments_tumor/modeled_segments
        - model_segments_tumor/modeled_segments_begin
        - model_segments_tumor/normal_het_allelic_counts
        - call_copy_ratio_segments_tumor/called_legacy_segments
        - call_copy_ratio_segments_tumor/called_segments
        - plot_denoised_copy_ratios_tumor/denoised_plot
        - plot_denoised_copy_ratios_tumor/delta_mad
        - plot_denoised_copy_ratios_tumor/denoised_mad
        - plot_denoised_copy_ratios_tumor/scaled_delta_mad
        - plot_denoised_copy_ratios_tumor/standardized_mad
        - plot_modeled_segments_tumor/output
    out: [output]

  tar_outputs_normal:
    run: ../tools/tar.cwl
    in:
      output_file:
        source: model_segments_normal/modeled_segments
        valueFrom: |
          $(self ? self.nameroot.split('.')[0]+".normal.gatk_cnv.tar.gz" : null)
      input_files:
        source:
        - collect_read_counts_normal/output
        - collect_allelic_counts_normal/output
        - denoise_read_counts_normal/denoised_copy_ratios
        - denoise_read_counts_normal/standardized_copy_ratios
        - model_segments_normal/allele_fraction_legacy_segments
        - model_segments_normal/allele_fraction_parameters
        - model_segments_normal/allele_fraction_parameters_begin
        - model_segments_normal/copy_ratio_legacy_segments
        - model_segments_normal/copy_ratio_only_segments
        - model_segments_normal/copy_ratio_parameters
        - model_segments_normal/copy_ratio_parameters_begin
        - model_segments_normal/het_allelic_counts
        - model_segments_normal/modeled_segments
        - model_segments_normal/modeled_segments_begin
        - model_segments_normal/normal_het_allelic_counts
        - call_copy_ratio_segments_normal/called_legacy_segments
        - call_copy_ratio_segments_normal/called_segments
        - plot_denoised_copy_ratios_normal/denoised_plot
        - plot_denoised_copy_ratios_normal/delta_mad
        - plot_denoised_copy_ratios_normal/denoised_mad
        - plot_denoised_copy_ratios_normal/scaled_delta_mad
        - plot_denoised_copy_ratios_normal/standardized_mad
        - plot_modeled_segments_normal/output
    out: [output]

$namespaces:
  sbg: https://sevenbridges.com
hints:
- class: 'sbg:maxNumberOfParallelInstances'
  value: 2
'sbg:license': Apache License 2.0
'sbg:publisher': KFDRC
'sbg:categories':
- BAM
- CNV
- COUNTS
- CRAM
- FUNCOTATOR
- HDF5
- NORMAL
- PON
- SEG
- SOMATIC
- TUMOR
'sbg:links':
- id: 'https://github.com/kids-first/kf-somatic-workflow/releases/tag/v0.2.0'
  label: github-release
