cwlVersion: v1.1
class: Workflow
id: kfdrc_gatk_create_cnv_pon_wf
label: KFDRC GATK Create CNV Panel of Normals Workflow
doc: |+
  # KFDRC GATK Create CNV Panel Of Normals Workflow

  ![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

  The Kids First Data Resource Center GATK Create CNV Panel of Normals Workflow is
  a direct port of the GATK best practices workflow WDL. Use this workflow for
  creating a GATK CNV Panel of Normals given a list of normal samples aligned
  reads. Supports both WGS and WES/WXS.

  ## Panel of Normal Recommendations

  ### Best Case Scenario

  Our Panel of Normal standards are the same as those defined in the Broad/GATK
  Documentation here:
  https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON-.

  In the best case scenario, you will have 40 samples derived from healthy
  tissue. Ideally these samples come from young and healthy individuals (thus
  eliminating the chance of a sample containing undiagnosed tumor tissue). Those
  40 samples will then be processed in the same way as the tumor sample.
  Processing includes all technical variables through sequencing (library
  preparation, sequencing platform, etc.). For WGS or WXS that spans the
  allosomes, the samples used for the panel of normals should have uniform sex.

  ### Imperfect Scenarios

  Often when working on a project, there will not be 40 normal samples available
  to the user. If multiple sequencing approaches/centers are used or you have
  allosomal regions of interest, the number of required normals can rapidly
  expand. If you are unable to acquire the 40 samples, it is still possible to
  create a smaller panel of normals. According to GATK, there is no definitive
  rule for the number of samples and even a small panel is better than no panel.

  Internally, we have observed that as the panel becomes smaller, more calls are
  made and those calls are less accurate. For WGS, we have observed that once
  panels become smaller than 15 samples, the calls begin to noticeably
  deteriorate. For WXS, calling deteriorates somewhere under 25 samples.

  ### Run at Own Risk

  It is possible to create a panel of normals from a single normal sample.
  Internally, we have observed that the calls generated from this scenario are
  highly erratic. In the case of WXS, the calls made had little overlap with the
  calls made from a larger panel. In the case of WGS, we had instances where the
  results looked slightly worse than the 15 sample panel of normals but we also
  had instances where the calls in no way resembled those from the larger panel.
  Again these results were all better than no panel but were far from what we
  would call reliable calls.

  #### Potential Workflow

  1. Your tumor sample shares an sequencing approach with 40 normals samples (same sex if the approach includes allosomal regions) in the same project.
  1. If not, reach out to the sequencing center that provided the tumor sample and request a set of normal samples that share the sequencing approach.
  1. If the center cannot provide such samples, search for public or similarly-controlled samples that share the sequencing approach. You will probably have most luck with WGS.

  At this point you have exhausted your sample sources. From here you can evaluate where you stand:

  - <span style="color:green;">`40 or more`</span>: you meet the GATK recommended minimum. Feel free to proceed.
  - <span style="color:yellowgreen;">`30 to 40`</span>: the results do not meet the GATK minimum but the panel is still rather sizable and worth a run
  - <span style="color:orange;">`20 to 30`</span>: calls from WXS or low-coverage WGS begin to deteriorate; run at your own risk!
  - <span style="color:orange;">`10 to 20`</span>: calls from high-coverage WGS begin to deteriorate; run at your own risk!
  - <span style="color:red">`1 to 10`</span>: do not trust calls without extensive follow-up analysis!!!


  ## Generalized steps:

  - prepares a genomic intervals list with PreprocessIntervals
  - collects read coverage counts across the preprocessed intervals
  - annotate the preprocessed intervals
  - creates a CNV PoN with CreateReadCountPanelOfNormals using read coverage
    counts

  ## Notes from GATK:

  - The intervals argument (`input_intervals` and/or `input_interval_list`) is
    required for both WGS and WES workflows and accepts formats compatible with
  the GATK -L argument. These intervals will be padded on both sides by the amount
  specified by `padding` (default 250) and split into bins of length specified by
  `bin_length` (default 1000; specify 0 to skip binning, e.g., for WES).  For WGS,
  the intervals should simply cover the autosomal chromosomes (sex chromosomes may
  be included, but care should be taken to 1) avoid creating panels of mixed sex,
  and 2) denoise case samples only with panels containing only individuals of the
  same sex as the case samples).
  - Intervals can be blacklisted from coverage collection and all downstream steps
    by using the blacklist intervals argument (`input_exclude_intervals` and/or
  `input_exclude_intervals_list`), which accepts formats compatible with the GATK
  -XL argument. This may be useful for excluding centromeric regions, etc. from
  analysis.  Alternatively, these regions may be manually filtered from the final
  callset.

  ### Other Resources

  dockerfiles: https://github.com/d3b-center/bixtools

requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: SubworkflowFeatureRequirement
inputs:
  input_aligned_reads: {type: 'File[]', secondaryFiles: [{pattern: ".bai", required: false},
      {pattern: "^.bai", required: false}, {pattern: ".crai", required: false}, {
        pattern: "^.crai", required: false}], doc: "Files containing aligned reads\
      \ from normal. Also include associated index.", 'sbg:fileTypes': "CRAM,BAM"}
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
  mappability_track: {type: 'File?', secondaryFiles: [{pattern: ".tbi", required: false}],
    doc: "Umap single-read mappability track. This is a BED file in .bed or .bed.gz\
      \ format that identifies uniquely mappable regions of the genome. The track\
      \ should correspond to the appropriate read length and overlapping intervals\
      \ must be merged. See https://bismap.hoffmanlab.org/.", 'sbg:fileTypes': "BED,BED.GZ"}
  segmental_duplication_track: {type: 'File?', secondaryFiles: [{pattern: ".tbi",
        required: false}], doc: "Segmental-duplication track. This is a BED file in\
      \ .bed or .bed.gz format that identifies segmental-duplication regions of the\
      \ genome. Overlapping intervals must be merged.", 'sbg:fileTypes': "BED,BED.GZ"}
  do_explicit_gc_correction: {type: 'boolean', doc: "If true, perform explicit GC-bias\
      \ correction when creating PoN and in subsequent denoising of case samples.\
      \ If false, rely on PCA-based denoising to correct for GC bias."}
  feature_query_lookahead: {type: 'int?', doc: "Number of bases to cache when querying\
      \ feature tracks."}
  readcount_format: {type: ['null', {type: enum, name: 'readcount_format', symbols: [
          "HDF5", "TSV"]}], doc: "Desired format of outputs from CollectReadCounts."}
  output_basename: {type: 'string?', doc: "String to use as basename for PoN and preprocessed\
      \ intervals files."}
  do_impute_zeros: {type: 'boolean?', doc: "If true, impute zero-coverage values as\
      \ the median of the non-zero values in the corresponding interval. (This is\
      \ applied after all filters.)"}
  extreme_sample_median_percentile: {type: 'float?', doc: "Samples with a median (across\
      \ genomic intervals) of fractional coverage normalized by genomic-interval medians\
      \ below this percentile or above the complementary percentile are filtered out.\
      \ (This is the fourth filter applied.)"}
  extreme_outlier_truncation_percentile: {type: 'float?', doc: "Fractional coverages\
      \ normalized by genomic-interval medians that are below this percentile or above\
      \ the complementary percentile are set to the corresponding percentile value.\
      \ (This is applied after all filters and imputation.)"}
  maximum_zeros_in_interval_percentage: {type: 'float?', doc: "Genomic intervals with\
      \ a fraction of zero-coverage samples above this percentage are filtered out.\
      \ (This is the third filter applied.)"}
  maximum_zeros_in_sample_percentage: {type: 'float?', doc: "Samples with a fraction\
      \ of zero-coverage genomic intervals above this percentage are filtered out.\
      \ (This is the second filter applied.)"}
  minimum_interval_median_percentile: {type: 'float?', doc: "Genomic intervals with\
      \ a median (across samples) of fractional coverage (optionally corrected for\
      \ GC bias) less than or equal to this percentile are filtered out. (This is\
      \ the first filter applied.)"}
  number_of_eigensamples: {type: 'int?', doc: "Number of eigensamples to use for truncated\
      \ SVD and to store in the panel of normals. The number of samples retained after\
      \ filtering will be used instead if it is smaller than this."}
  create_pon_max_memory: {type: 'int?', doc: "Max memory in GB to allocate to the\
      \ task"}
  create_pon_cpus: {type: 'int?', doc: "Cores to allocate to the task"}
outputs:
  preprocessed_intervals: {type: File, outputSource: preprocess_intervals/output,
    doc: "Preprocessed Picard interval-list file used to collect read counts for PoN."}
  read_counts: {type: 'File[]', outputSource: collect_read_counts/output, doc: "Counts\
      \ files generated for each input aligned reads file. Used to generate PoN."}
  cnv_pon: {type: File, outputSource: create_read_count_panel_of_normals/output, doc: "Panel-of-normals\
      \ file in HDF5 format."}
steps:
  preprocess_intervals:
    run: ../tools/gatk_preprocessintervals.cwl
    in:
      reference: reference_fasta
      sequence_dictionary: reference_dict
      input_intervals: input_intervals
      input_interval_list: input_interval_list
      input_exclude_intervals: input_exclude_intervals
      input_exclude_interval_list: input_exclude_interval_list
      padding: padding
      bin_length: bin_length
      interval_merging_rule: {default: "OVERLAPPING_ONLY"}
      output_prefix:
        source: output_basename
        valueFrom: $(self).gatk_cnv
    out: [output]
  annotate_intervals:
    run: ../tools/gatk_annotateintervals.cwl
    in:
      do_explicit_gc_correction: do_explicit_gc_correction
      reference: reference_fasta
      sequence_dictionary: reference_dict
      input_interval_list: preprocess_intervals/output
      mappability_track: mappability_track
      segmental_duplication_track: segmental_duplication_track
      interval_merging_rule: {default: "OVERLAPPING_ONLY"}
      feature_query_lookahead: feature_query_lookahead
    out: [output]
  collect_read_counts:
    run: ../tools/gatk_collectreadcounts.cwl
    hints:
    - class: sbg:AWSInstanceType
      value: c5.9xlarge
    scatter: input_aligned_reads
    in:
      reference: reference_fasta
      sequence_dictionary: reference_dict
      input_interval_list: preprocess_intervals/output
      input_aligned_reads: input_aligned_reads
      interval_merging_rule: {default: "OVERLAPPING_ONLY"}
      output_format: readcount_format
    out: [output]
  create_read_count_panel_of_normals:
    run: ../tools/gatk_createreadcountpanelofnormals.cwl
    in:
      input_counts: collect_read_counts/output
      input_annotated_intervals: annotate_intervals/output
      output_basename:
        source: output_basename
        valueFrom: $(self).gatk_cnv
      do_impute_zeros: do_impute_zeros
      extreme_sample_median_percentile: extreme_sample_median_percentile
      extreme_outlier_truncation_percentile: extreme_outlier_truncation_percentile
      maximum_zeros_in_interval_percentage: maximum_zeros_in_interval_percentage
      maximum_zeros_in_sample_percentage: maximum_zeros_in_sample_percentage
      minimum_interval_median_percentile: minimum_interval_median_percentile
      number_of_eigensamples: number_of_eigensamples
      max_memory: create_pon_max_memory
      cpus: create_pon_cpus
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
- HDF5
- PON
'sbg:links':
- id: 'https://github.com/kids-first/kf-somatic-workflow/releases/tag/v4.0.0'
  label: github-release
