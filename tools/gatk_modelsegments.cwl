cwlVersion: v1.0
class: CommandLineTool
id: gatk_modelsegments
doc: "Models segmented copy ratios from denoised read counts and segmented minor-allele fractions from allelic counts"
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.max_memory*1000)
  - class: DockerRequirement
    dockerPull: 'broadinstitute/gatk:4.2.4.1'
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      $( inputs.denoised_copy_ratios ? 'gatk' : 'echo gatk' )
  - position: 1
    shellQuote: false
    prefix: "--java-options"
    valueFrom: >-
      $("\"-Xmx"+Math.floor(inputs.max_memory*1000/1.074 - 1)+"M\"")
  - position: 2
    shellQuote: false
    valueFrom: ModelSegments
  - position: 3
    shellQuote: false
    prefix: "--output"
    valueFrom: >-
      .
  - position: 3
    shellQuote: false
    prefix: "--output-prefix"
    valueFrom: >-
      ${var pre = inputs.output_prefix ? inputs.output_prefix : inputs.denoised_copy_ratios ? inputs.denoised_copy_ratios.nameroot.split('.').slice(0,-1).join('.') : inputs.allelic_counts ? inputs.allelic_counts.nameroot.split('.').slice(0,-1).join('.') : 'output_ModelSegments'; return pre}
inputs:
  allelic_counts:
    type: 'File?'
    doc: "Input file containing allelic counts (output of CollectAllelicCounts)."
    inputBinding:
      position: 3
      prefix: "--allelic-counts"
  denoised_copy_ratios:
    type: 'File?'
    doc: "Input file containing denoised copy ratios (output of DenoiseReadCounts)."
    inputBinding:
      position: 3
      prefix: "--denoised-copy-ratios"
  genotyping_base_error_rate:
    type: 'float?'
    doc: "Maximum base-error rate for genotyping and filtering homozygous allelic counts, if available. The likelihood for an allelic count to be generated from a homozygous site will be integrated from zero base-error rate up to this value. Decreasing this value will increase the number of sites assumed to be heterozygous for modeling."
    inputBinding:
      position: 3
      prefix: "--genotyping-base-error-rate"
  genotyping_homozygous_log_ratio_threshold:
    type: 'float?'
    doc: "Log-ratio threshold for genotyping and filtering homozygous allelic counts, if available. Increasing this value will increase the number of sites assumed to be heterozygous for modeling."
    inputBinding:
      position: 3
      prefix: "--genotyping-homozygous-log-ratio-threshold"
  kernel_approximation_dimension:
    type: 'int?'
    doc: "Dimension of the kernel approximation. A subsample containing this number of data points will be used to construct the approximation for each chromosome. If the total number of data points in a chromosome is greater than this number, then all data points in the chromosome will be used. Time complexity scales quadratically and space complexity scales linearly with this parameter."
    inputBinding:
      position: 3
      prefix: "--kernel-approximation-dimension"
  kernel_scaling_allele_fraction:
    type: 'float?'
    doc: "Relative scaling S of the kernel K_AF for allele-fraction segmentation to the kernel K_CR for copy-ratio segmentation. If multidimensional segmentation is performed, the total kernel used will be K_CR + S * K_AF."
    inputBinding:
      position: 3
      prefix: "--kernel-scaling-allele-fraction"
  kernel_variance_allele_fraction:
    type: 'float?'
    doc: "Variance of Gaussian kernel for allele-fraction segmentation, if performed. If zero, a linear kernel will be used."
    inputBinding:
      position: 3
      prefix: "--kernel-variance-allele-fraction"
  kernel_variance_copy_ratio:
    type: 'float?'
    doc: "Variance of Gaussian kernel for copy-ratio segmentation, if performed. If zero, a linear kernel will be used."
    inputBinding:
      position: 3
      prefix: "--kernel-variance-copy-ratio"
  maximum_number_of_segments_per_chromosome:
    type: 'int?'
    doc: "Maximum number of segments allowed per chromosome."
    inputBinding:
      position: 3
      prefix: "--maximum-number-of-segments-per-chromosome"
  maximum_number_of_smoothing_iterations:
    type: 'int?'
    doc: "Maximum number of iterations allowed for segmentation smoothing."
    inputBinding:
      position: 3
      prefix: "--maximum-number-of-smoothing-iterations"
  minimum_total_allele_count_case:
    type: 'int?'
    doc: "Minimum total count for filtering allelic counts in the case sample, if available. The default value of zero is appropriate for matched-normal mode; increase to an appropriate value for case-only mode."
    inputBinding:
      position: 3
      prefix: "--minimum-total-allele-count-case"
  minimum_total_allele_count_normal:
    type: 'int?'
    doc: "Minimum total count for filtering allelic counts in the matched-normal sample, if available."
    inputBinding:
      position: 3
      prefix: "--minimum-total-allele-count-normal"
  minor_allele_fraction_prior_alpha:
    type: 'float?'
    doc: "Alpha hyperparameter for the 4-parameter beta-distribution prior on segment minor-allele fraction. The prior for the minor-allele fraction f in each segment is assumed to be Beta(alpha, 1, 0, 1/2). Increasing this hyperparameter will reduce the effect of reference bias at the expense of sensitivity."
    inputBinding:
      position: 3
      prefix: "--minor-allele-fraction-prior-alpha"
  normal_allelic_counts:
    type: 'File?'
    doc: "Input file containing allelic counts for a matched normal (output of CollectAllelicCounts)."
    inputBinding:
      position: 3
      prefix: "--normal-allelic-counts"
  number_of_burn_in_samples_allele_fraction:
    type: 'int?'
    doc: "Number of burn-in samples to discard for allele-fraction model."
    inputBinding:
      position: 3
      prefix: "--number-of-burn-in-samples-allele-fraction"
  number_of_burn_in_samples_copy_ratio:
    type: 'int?'
    doc: "Number of burn-in samples to discard for copy-ratio model."
    inputBinding:
      position: 3
      prefix: "--number-of-burn-in-samples-copy-ratio"
  number_of_changepoints_penalty_factor:
    type: 'int?'
    doc: "Factor A for the penalty on the number of changepoints per chromosome for segmentation. Adds a penalty of the form A * C * [1 + log (N / C)], where C is the number of changepoints in the chromosome, to the cost function for each chromosome. Must be non-negative."
    inputBinding:
      position: 3
      prefix: "--number-of-changepoints-penalty-factor"
  number_of_samples_allele_fraction:
    type: 'int?'
    doc: "Total number of MCMC samples for allele-fraction model."
    inputBinding:
      position: 3
      prefix: "--number-of-samples-allele-fraction"
  number_of_samples_copy_ratio:
    type: 'int?'
    doc: "Total number of MCMC samples for copy-ratio model."
    inputBinding:
      position: 3
      prefix: "--number-of-samples-copy-ratio"
  number_of_smoothing_iterations_per_fit:
    type: 'int?'
    doc: "Number of segmentation-smoothing iterations per MCMC model refit. (Increasing this will decrease runtime, but the final number of segments may be higher. Setting this to 0 will completely disable model refitting between iterations.)"
    inputBinding:
      position: 3
      prefix: "--number-of-smoothing-iterations-per-fit"
  smoothing_credible_interval_threshold_allele_fraction:
    type: 'float?'
    doc: "Number of 10% equal-tailed credible-interval widths to use for allele-fraction segmentation smoothing."
    inputBinding:
      position: 3
      prefix: "--smoothing-credible-interval-threshold-allele-fraction"
  smoothing_credible_interval_threshold_copy_ratio:
    type: 'float?'
    doc: "Number of 10% equal-tailed credible-interval widths to use for copy-ratio segmentation smoothing."
    inputBinding:
      position: 3
      prefix: "--smoothing-credible-interval-threshold-copy-ratio"
  window_size:
    type:
      - 'null'
      - type: array
        items: int
        inputBinding:
          prefix: "--window-size"
    doc: "Window sizes to use for calculating local changepoint costs. For each window size, the cost for each data point to be a changepoint will be calculated assuming that the point demarcates two adjacent segments of that size. Including small (large) window sizes will increase sensitivity to small (large) events. Duplicate values will be ignored."
    inputBinding:
      position: 3
  output_prefix:
    type: 'string?'
    doc: "String to use as the prefix for the outputs."
  max_memory:
    type: 'int?'
    default: 20
    doc: "Maximum GB of RAM to allocate for this tool."
outputs:
  allele_fraction_legacy_segments: { type: 'File?', outputBinding: { glob: "*.af.igv.seg" } }
  allele_fraction_parameters: { type: 'File?', outputBinding: { glob: "*.modelFinal.af.param" } }
  allele_fraction_parameters_begin: { type: 'File?', outputBinding: { glob: "*.modelBegin.af.param" } }
  copy_ratio_legacy_segments: { type: 'File?', outputBinding: { glob: "*.cr.igv.seg" } }
  copy_ratio_only_segments: { type: 'File?', outputBinding: { glob: "*.cr.seg" } }
  copy_ratio_parameters: { type: 'File?', outputBinding: { glob: "*.modelFinal.cr.param" } }
  copy_ratio_parameters_begin: { type: 'File?', outputBinding: { glob: "*.modelBegin.cr.param" } }
  het_allelic_counts: { type: 'File?', outputBinding: { glob: "*.hets.tsv" } }
  modeled_segments: { type: 'File?', outputBinding: { glob: "*.modelFinal.seg" } }
  modeled_segments_begin: { type: 'File?', outputBinding: { glob: "*.modelBegin.seg" } }
  normal_het_allelic_counts: { type: 'File?', outputBinding: { glob: "*.hets.normal.tsv" } }
