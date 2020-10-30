cwlVersion: v1.0
class: CommandLineTool
id: canvas
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/canvas:1.11.0'
  - class: ResourceRequirement
    ramMin: 32000
    coresMin: 16
  - class: InlineJavascriptRequirement
arguments: 
  - position: 1
    shellQuote: false
    valueFrom: >-
      ln -s $(inputs.control_bam.path) .; ln -s $(inputs.control_bam.secondaryFiles[0].path) ./$(inputs.control_bam.basename).bai

      ln -s $(inputs.tumor_bam.path) .; ln -s $(inputs.tumor_bam.secondaryFiles[0].path) ./$(inputs.tumor_bam.basename).bai

      mono /1.11.0/Canvas.exe Somatic-Enrichment
      -b $(inputs.tumor_bam.basename)
      --manifest=$(inputs.manifest.path)
      --control-bam=$(inputs.control_bam.basename)
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
  tumor_bam: {type: File, doc: "tumor bam file", secondaryFiles: [.bai]}
  manifest: {type: File, doc: "Nextera manifest file"}
  control_bam: {type: File, doc: "Bam file of unmatched control sample", secondaryFiles: [.bai]}
  b_allele_vcf: {type: File, doc: "vcf containing SNV b-alleles sites (only sites with PASS will be used)"}
  sample_name: string
  reference: {type: File, doc: "Canvas-ready kmer file"}
  genomeSize_file: {type: File, doc: "GenomeSize.xml"}
  genome_fasta: {type: File, doc: "Genome.fa", secondaryFiles: [.fai]}
  filter_bed: {type: File, doc: "bed file of regions to skip"}
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

  
