cwlVersion: v1.0
class: CommandLineTool
id: consensus_merge
doc: >-
  This tool will run a Python script that merges somatic calls from several callers
      into a single consensus VCF
requirements:
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/consensus-merge:1.1.0'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: ${ return inputs.ram * 1000 }
    coresMin: $(inputs.cores)
  - class: ShellCommandRequirement

baseCommand: [/usr/bin/consensus_merge.py]

inputs:
  strelka2_vcf:
    type: File
    inputBinding:
      position: 1
      prefix: '--strelka2_vcf'
    secondaryFiles: ['.tbi']
    doc: 'Strelka2 VCF with MNPs fixed'
  mutect2_vcf:
    type: File
    inputBinding:
      position: 2
      prefix: '--mutect2_vcf'
    secondaryFiles: ['.tbi']
    doc: 'Mutect2 VCF'
  lancet_vcf:
    type: File
    inputBinding:
      position: 3
      prefix: '--lancet_vcf'
    secondaryFiles: ['.tbi']
    doc: 'Lancet VCF'
  vardict_vcf:
    type: File
    inputBinding:
      position: 4
      prefix: '--vardict_vcf'
    secondaryFiles: ['.tbi']
    doc: 'VarDict VCF'
  cram:
    type: File
    inputBinding:
      position: 5
      prefix: '--cram'
    secondaryFiles: ['.crai']
    doc: 'CRAM or BAM file'
  reference:
    type: File
    inputBinding:
      position: 6
      prefix: '--reference'
    secondaryFiles: ['.fai']
    doc: 'Path to FASTA to which CRAM is aligned'
  ncallers:
    type: int?
    inputBinding:
      position: 7
      prefix: '--ncallers'
    doc: 'Optional number of callers required for consensus'
  output_basename:
    type: string
    inputBinding:
      position: 8
      prefix: '--output_basename'
  hotspot_source:
    type: string?
    inputBinding:
      position: 9
      prefix: '--hotspot_source'
  contig_bed:
    type: File?
    inputBinding:
      position: 10
      prefix: '--contig_bed'
  cores:
    type: int?
    default: 16
  ram:
    type: int?
    default: 3 
    doc: 'RAM requirement in GB'

outputs:
  output:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: ['.tbi']
