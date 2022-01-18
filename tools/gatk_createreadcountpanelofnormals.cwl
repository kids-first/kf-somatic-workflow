cwlVersion: v1.0
class: CommandLineTool
id: gatk_createreadcountpanelofnormals
doc: "Creates a panel of normals (PoN) for read-count denoising given the read counts for samples in the panel"
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.max_memory * 1000)
    coresMin: $(inputs.cpus)
  - class: DockerRequirement
    dockerPull: 'broadinstitute/gatk:4.2.4.1'
baseCommand: [gatk]
arguments:
  - position: 0
    shellQuote: false
    prefix: "--java-options"
    valueFrom: >-
      '-Xms$(inputs.max_memory*1000/4)M -Xmx$(Math.floor(inputs.max_memory*1000/1.074 - 1))M'
  - position: 1
    shellQuote: false
    valueFrom: >-
      CreateReadCountPanelOfNormals 
  - position: 2
    shellQuote: false
    prefix: "-O"
    valueFrom: >-
      ${if(inputs.output_basename) {return inputs.output_basename+"."} else {return "" } }cnv.pon.hdf5
inputs:
  input_counts:
    type:
      type: array
      items: File
      inputBinding:
        prefix: "-I"
    doc: "Input TSV or HDF5 files containing integer read counts in genomic intervals for all samples in the panel of normals (output of CollectReadCounts). Intervals must be identical and in the same order for all samples."
    inputBinding:
      position: 2
  input_annotated_intervals:
    type: 'File?'
    doc: "Input file containing annotations for GC content in genomic intervals (output of AnnotateIntervals). If provided, explicit GC correction will be performed before performing SVD. Intervals must be identical to and in the same order as those in the input read-counts files."
    inputBinding:
      position: 2
      prefix: "--annotated-intervals"
  output_basename:
    type: 'string?'
    doc: "String to use as the prefix for the outputs."
  do_impute_zeros:
    type: 'boolean?'
    doc: "If true, impute zero-coverage values as the median of the non-zero values in the corresponding interval. (This is applied after all filters.)"
    inputBinding:
      position: 2
      prefix: "--do-impute-zeros"
  extreme_outlier_truncation_percentile:
    type: 'float?'
    doc: "Fractional coverages normalized by genomic-interval medians that are below this percentile or above the complementary percentile are set to the corresponding percentile value. (This is applied after all filters and imputation.)"
    inputBinding:
      position: 2
      prefix: "--extreme-outlier-truncation-percentile"
  extreme_sample_median_percentile:
    type: 'float?'
    doc: "Samples with a median (across genomic intervals) of fractional coverage normalized by genomic-interval medians below this percentile or above the complementary percentile are filtered out. (This is the fourth filter applied.)"
    inputBinding:
      position: 2
      prefix: "--extreme-sample-median-percentile"
  maximum_zeros_in_interval_percentage:
    type: 'float?'
    doc: "Genomic intervals with a fraction of zero-coverage samples above this percentage are filtered out. (This is the third filter applied.)"
    inputBinding:
      position: 2
      prefix: "--maximum-zeros-in-interval-percentage"
  maximum_zeros_in_sample_percentage:
    type: 'float?'
    doc: "Samples with a fraction of zero-coverage genomic intervals above this percentage are filtered out. (This is the second filter applied.)"
    inputBinding:
      position: 2
      prefix: "--maximum-zeros-in-sample-percentage"
  minimum_interval_median_percentile:
    type: 'float?'
    doc: "Genomic intervals with a median (across samples) of fractional coverage (optionally corrected for GC bias) less than or equal to this percentile are filtered out. (This is the first filter applied.)"
    inputBinding:
      position: 2
      prefix: "--minimum-interval-median-percentile"
  number_of_eigensamples:
    type: 'int?'
    doc: "Number of eigensamples to use for truncated SVD and to store in the panel of normals. The number of samples retained after filtering will be used instead if it is smaller than this."
    inputBinding:
      position: 2
      prefix: "--number-of-eigensamples"
  max_memory: { type: 'int?', default: 16, doc: "Max memory in GB to allocate to the task" }
  cpus: { type: 'int?', default: 8, doc: "Cores to allocate to the task" }
outputs:
  output: { type: File, outputBinding: { glob: "*.hdf5" } }
