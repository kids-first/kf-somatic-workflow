cwlVersion: v1.0
class: CommandLineTool
id: canvas
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/canvas:1.11.0'
  - class: ResourceRequirement
    ramMin: 32000
    coresMin: 16
  - class: InlineJavascriptRequirement
baseCommand: [mono]
arguments: 
  - position: 1
    shellQuote: false
    valueFrom: >-
      /1.11.0/Canvas.exe Somatic-Enrichment
      -b $(inputs.tumor_bam.path)
      --manifest=$(inputs.manifest.path)
      --control-bam=$(inputs.control_bam.path)
      --b-allele-vcf=$(inputs.b_allele_vcf.path)
      --exclude-non-het-b-allele-sites
      --sample-name=$(inputs.sample_name)
      --genome-folder=$(inputs.genome_fasta.dirname)
      -o ./
      -r $(inputs.reference.path)
      --filter-bed=$(inputs.filter_bed.path)

      mv CNV.vcf.gz $(inputs.output_basename).canvas.CNV.vcf.gz &&
      tabix $(inputs.output_basename).canvas.CNV.vcf.gz

      mv CNV.CoverageAndVariantFrequency.txt $(inputs.output_basename).canvas.CNV.CoverageAndVariantFrequency.txt

      tar -czf TempCNV_$(inputs.sample_name).tar.gz
      TempCNV_$(inputs.sample_name)

inputs:
  tumor_bam: {type: File, label: tumor bam file, secondaryFiles: [.bai]}
  manifest: {type: File, label: Nextera manifest file}
  control_bam: {type: ['null', File], label: Bam file of unmatched control sample (optional), secondaryFiles: [.bai]}
  b_allele_vcf: {type: File, label: vcf containing SNV b-alleles sites (only sites with PASS will be used)}
  sample_name: string
  reference: {type: File, label: Canvas-ready kmer file}
  genomeSize_file: {type: File, label: GenomeSize.xml}
  genome_fasta: {type: File, label: Genome.fa, secondaryFiles: [.fai]}
  filter_bed: {type: File, label: bed file of regions to skip }
  output_basename: string

outputs:
  output_vcf:
    type: File
    outputBinding:
      glob: '*.CNV.vcf.gz'
    secondaryFiles: [.tbi]
  output_txt:
    type: File
    outputBinding:
      glob: '*.CNV.CoverageAndVariantFrequency.txt'
  output_folder:
    type: File
    outputBinding:
      glob: '*.tar.gz'

  