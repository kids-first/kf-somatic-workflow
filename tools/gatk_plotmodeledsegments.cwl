cwlVersion: v1.0
class: CommandLineTool
id: gatk_plotmodelsegments
doc: "Creates plots of denoised and segmented copy-ratio and minor-allele-fraction estimates"
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
    valueFrom: PlotModeledSegments
  - position: 2
    shellQuote: false
    prefix: "--output"
    valueFrom: >-
      .
  - position: 2
    shellQuote: false
    prefix: "--output-prefix"
    valueFrom: >-
      ${var pre = inputs.output_prefix ? inputs.output_prefix : inputs.segments ? inputs.segments.nameroot : 'plotsegments'; return pre}
inputs:
  allelic_counts:
    type: 'File?'
    doc: "Input file containing allelic counts at heterozygous sites (.hets.tsv output of ModelSegments)."
    inputBinding:
      position: 2
      prefix: "--allelic-counts"
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
  segments:
    type: 'File?'
    doc: "Input file containing modeled segments (output of ModelSegments)."
    inputBinding:
      position: 2
      prefix: "--segments"
  sequence_dictionary:
    type: 'File'
    doc: "Use the given sequence dictionary as the master/canonical sequence dictionary. Must be a .dict file."
    inputBinding:
      position: 2
      prefix: "--sequence-dictionary"
  output_prefix:
    type: 'string?'
    doc: "String to use as the prefix for the outputs."
outputs:
  output: { type: 'File?', outputBinding: { glob: "*.modeled.png" } }
