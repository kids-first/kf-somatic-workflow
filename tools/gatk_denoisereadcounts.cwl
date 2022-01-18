cwlVersion: v1.0
class: CommandLineTool
id: gatk_denoisereadcounts
doc: >-
  Denoises read counts to produce denoised copy ratios.
  
  Typically, a panel of normals produced by CreateReadCountPanelOfNormals is provided as input. The
  input counts are then standardized by 1) transforming to fractional coverage, 2) performing
  optional explicit GC-bias correction (if the panel contains GC-content annotated intervals), 3)
  filtering intervals to those contained in the panel, 4) dividing by interval medians contained in
  the panel, 5) dividing by the sample median, and 6) transforming to log2 copy ratio. The result
  is then denoised by subtracting the projection onto the specified number of principal components
  from the panel.
  
  If no panel is provided, then the input counts are instead standardized by 1) transforming to
  fractional coverage, 2) performing optional explicit GC-bias correction (if GC-content annotated
  intervals are provided), 3) dividing by the sample median, and 4) transforming to log2 copy ratio.
  No denoising is performed, so the denoised result is simply taken to be identical to the
  standardized result.
  
  If performed, explicit GC-bias correction is done by GCBiasCorrector.
  
  Note that number-of-eigensamples principal components from the input panel will be used for
  denoising; if only fewer are available in the panel, then they will all be used. This parameter
  can thus be used to control the amount of denoising, which will ultimately affect the sensitivity
  of the analysis.
  
  See comments for CreateReadCountPanelOfNormals regarding coverage on sex chromosomes. If sex
  chromosomes are not excluded from coverage collection, it is strongly recommended that case
  samples are denoised only with panels containing only individuals of the same sex as the case
  samples. 
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'broadinstitute/gatk:4.2.4.1'
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      $( inputs.read_counts ? 'gatk' : 'echo gatk' )
  - position: 1
    shellQuote: false
    valueFrom: DenoiseReadCounts 
  - position: 2
    shellQuote: false
    prefix: "--denoised-copy-ratios"
    valueFrom: >-
      ${var pre = inputs.output_prefix ? inputs.output_prefix : inputs.read_counts ? inputs.read_counts.nameroot : 'output'; return pre+'.denoisedCR.tsv'}
  - position: 2
    shellQuote: false
    prefix: "--standardized-copy-ratios"
    valueFrom: >-
      ${var pre = inputs.output_prefix ? inputs.output_prefix : inputs.read_counts ? inputs.read_counts.nameroot : 'output'; return pre+'.standardizedCR.tsv'}
inputs:
  annotated_intervals:
    type: 'File?'
    doc: "Input file containing annotations for GC content in genomic intervals (output of AnnotateIntervals). Intervals must be identical to and in the same order as those in the input read-counts file. If a panel of normals is provided, this input will be ignored."
    inputBinding:
      position: 2
      prefix: "--annotated-intervals"
  count_panel_of_normals:
    type: 'File?'
    doc: "Input HDF5 file containing the panel of normals (output of CreateReadCountPanelOfNormals)."
    inputBinding:
      position: 2
      prefix: "--count-panel-of-normals"
  read_counts:
    type: 'File?'
    doc: "Input TSV or HDF5 file containing integer read counts in genomic intervals for a single case sample (output of CollectReadCounts)."
    inputBinding:
      position: 2
      prefix: "-I"
  number_of_eigensamples:
    type: 'int?'
    doc: "Number of eigensamples to use for denoising. If not specified or if the number of eigensamples available in the panel of normals is smaller than this, all eigensamples will be used."
    inputBinding:
      position: 2
      prefix: "--number-of-eigensamples"
  output_prefix:
    type: 'string?'
    doc: "String to use as the prefix for the outputs."
outputs:
  denoised_copy_ratios: { type: 'File?', outputBinding: { glob: "*.denoisedCR.tsv" } }
  standardized_copy_ratios: { type: 'File?', outputBinding: { glob: "*.standardizedCR.tsv" } }
