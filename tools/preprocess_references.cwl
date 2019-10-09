cwlVersion: v1.0
class: CommandLineTool
id: preprocess_references
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/samtools:1.9'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 16000
    coresMin: 8
  
baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      cp $(inputs.input_fasta_file.path) .

      samtools index -@ 8 ./$(inputs.input_fasta_file.basename)

      cp $(inputs.af_only_gnomad_vcf.path) .

      tabix $(inputs.af_only_gnomad_vcf.basename)

      cp $(inputs.exac_common_vcf.path) .

      tabix $(inputs.exac_common_vcf.basename)

      cp $(inputs.hg38_strelka_bed.path) .

      tabix $(inputs.hg38_strelka_bed.basename)

inputs:
  input_fasta_file: File
  af_only_gnomad_vcf: File
  exac_common_vcf: File
  hg38_strelka_bed
outputs:
  indexed_reference_fasta:
    type: File
    outputBinding:
      glob: '$(inputs.input_fasta_file.basename)'
    secondaryFiles: [^.fai]
  indexed_af_only_gnomad_vcf:
    type: File
    outputBinding:
      glob: '$(inputs.af_only_gnomad_vcf.basename)'
    secondaryFiles: [.tbi]
  indexed_exac_common_vcf:
    type: File
    outputBinding:
      glob: '$(inputs.exac_common_vcf.basename)'
    secondaryFiles: [.tbi]
  indexed_hg38_strelka_bed:
    type: File
    outputBinding:
      glob: '$(inputs.hg38_strelka_bed.basename)'
    secondaryFiles: [.tbi]
