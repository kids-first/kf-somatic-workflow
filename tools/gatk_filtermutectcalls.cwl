cwlVersion: v1.0
class: CommandLineTool
id: gatk4_filtermutect2calls
label: GATK Filter Mutect2
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: ${ return inputs.max_memory * 1000 }
    coresMin: 2
baseCommand: [/gatk, FilterMutectCalls]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m"
      -V $(inputs.mutect_vcf.path)
      -O $(inputs.output_basename).mutect2_filtered.vcf.gz
      -R $(inputs.reference.path)
      --contamination-table $(inputs.contamination_table.path)
      --tumor-segmentation $(inputs.segmentation_table.path)
      --ob-priors $(inputs.ob_priors.path)
      --filtering-stats $(inputs.output_basename).mutect2_filtered.txt
      --stats $(inputs.mutect_stats.path)

inputs:
  mutect_vcf: {type: File, secondaryFiles: [.tbi]}
  mutect_stats: File
  reference: File
  output_basename: string
  contamination_table: File
  segmentation_table: File
  ob_priors: File
  max_memory: {type: int?, default: 4, doc: "GB of memory to allocate to the task"}

outputs:
  stats_table:
    type: File
    outputBinding:
      glob: '*.mutect2_filtered.txt'

  filtered_vcf:
    type: File
    outputBinding:
      glob: '*.mutect2_filtered.vcf.gz'
    secondaryFiles: ['.tbi']
