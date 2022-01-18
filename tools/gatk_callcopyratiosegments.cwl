cwlVersion: v1.0
class: CommandLineTool
id: gatk_callcopyratiosegments
doc: "Calls copy-ratio segments as amplified, deleted, or copy-number neutral."
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
      $( inputs.copy_ratio_segments ? 'gatk' : 'echo gatk' )
  - position: 1
    shellQuote: false
    valueFrom: CallCopyRatioSegments
  - position: 2
    shellQuote: false
    prefix: "-O"
    valueFrom: >-
      ${var pre = inputs.output_prefix ? inputs.output_prefix : inputs.copy_ratio_segments ? inputs.copy_ratio_segments.nameroot : 'output'; return pre+'.called.seg'}

inputs:
  calling_copy_ratio_z_score_threshold:
    type: 'float?'
    doc: "Threshold on z-score of non-log2 copy ratio used for calling segments."
    inputBinding:
      position: 2
      prefix: "--calling-copy-ratio-z-score-threshold"
  copy_ratio_segments:
    type: 'File?'
    doc: "Input file containing copy-ratio segments (.cr.seg output of ModelSegments)."
    inputBinding:
      position: 2
      prefix: "-I"
  neutral_segment_copy_ratio_lower_bound:
    type: 'float?'
    doc: "Lower bound on non-log2 copy ratio used for determining copy-neutral segments."
    inputBinding:
      position: 2
      prefix: "--neutral-segment-copy-ratio-lower-bound"
  neutral_segment_copy_ratio_upper_bound:
    type: 'float?'
    doc: "Upper bound on non-log2 copy ratio used for determining copy-neutral segments."
    inputBinding:
      position: 2
      prefix: "--neutral-segment-copy-ratio-upper-bound"
  outlier_neutral_segment_copy_ratio_z_score_threshold:
    type: 'float?'
    doc: "Threshold on z-score of non-log2 copy ratio used for determining outlier copy-neutral segments. If non-log2 copy ratio z-score is above this threshold for a copy-neutral segment, it is considered an outlier and not used in the calculation of the length-weighted mean and standard deviation used for calling."
    inputBinding:
      position: 2
      prefix: "--outlier-neutral-segment-copy-ratio-z-score-threshold"
  output_prefix:
    type: 'string?'
    doc: "String to use as the prefix for the outputs."

outputs:
  called_legacy_segments: { type: 'File?', outputBinding: { glob: "*called.igv.seg" } }
  called_segments: { type: 'File?', outputBinding: { glob: "*.called.seg" } }
