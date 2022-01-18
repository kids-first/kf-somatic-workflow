cwlVersion: v1.0
class: CommandLineTool
id: gatk_plotdenoisedcopyratios
doc: "Creates plots of denoised copy ratios"
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: ResourceRequirement
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
    valueFrom: PlotDenoisedCopyRatios 
  - position: 2
    shellQuote: false
    prefix: "--output"
    valueFrom: >-
      .
  - position: 2
    shellQuote: false
    prefix: "--output-prefix"
    valueFrom: >-
      ${var pre = inputs.output_prefix ? inputs.output_prefix : inputs.denoised_copy_ratios ? inputs.denoised_copy_ratios.nameroot.split('.').slice(0,-1).join('.') : 'plotdenoisecr'; return pre}
inputs:
  denoised_copy_ratios:
    type: 'File?'
    doc: "Input file containing denoised copy ratios (output of DenoiseReadCounts)."
    inputBinding:
      position: 2
      prefix: "--denoised-copy-ratios"
  minimum_contig_length:
    type: 'int?'
    doc: "Threshold length (in bp) for contigs to be plotted. Contigs with lengths less than this threshold will not be plotted. This can be used to filter out mitochondrial contigs, unlocalized contigs, etc."
    inputBinding:
      position: 2
      prefix: "--minimum-contig-length"
  sequence_dictionary:
    type: 'File'
    doc: "Use the given sequence dictionary as the master/canonical sequence dictionary. Must be a .dict file."
    inputBinding:
      position: 2
      prefix: "--sequence-dictionary"
  standardized_copy_ratios:
    type: 'File?'
    doc: "Input file containing standardized copy ratios (output of DenoiseReadCounts)."
    inputBinding:
      position: 2
      prefix: "--standardized-copy-ratios"
  output_prefix:
    type: 'string?'
    doc: "String to use as the prefix for the outputs."
outputs:
  denoised_plot: { type: 'File?', outputBinding: { glob: "*.denoised.png" } }
  delta_mad: { type: 'File?', outputBinding: { glob: "*.deltaMAD.txt" } }
  denoised_mad: { type: 'File?', outputBinding: { glob: "*.denoisedMAD.txt" } }
  scaled_delta_mad: { type: 'File?', outputBinding: { glob: "*.scaledDeltaMAD.txt" } }
  standardized_mad: { type: 'File?', outputBinding: { glob: "*.standardizedMAD.txt" } }
